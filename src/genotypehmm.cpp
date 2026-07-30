#include <stdexcept>
#include <cassert>
#include <limits>
#include <fstream>
#include <array>
#include <algorithm>
#include <cmath>
#include <vector>

#include "genotypehmm.h"
#include "bipartitioniterator.h"
#include "transitionprobabilitycomputer.h"
#include "binomial.h"

using namespace std;

GenotypeHMM::GenotypeHMM(
	ReadSet* read_set, 
	const uint32_t ploidy, 
	const float& recombrate,
	const float& eff_pop_size,
	const uint32_t& num_haplotypes, 
	vector<variant_information_t>* variant_info_table)
	:read_set(read_set),
	ploidy(ploidy),
	recombrate(recombrate),
	eff_pop_size(eff_pop_size),
	column_iterator(*read_set, variant_info_table),
	scaling_parameters(column_iterator.get_column_count(),-1.0L),
	num_haplotypes(num_haplotypes),
	variant_info_table(variant_info_table)
{
	//compute forward and backward probabilities
	std::cerr << "[Core::Genotyping] Initializing index structure for HMM." << std::endl;
	compute_index();
	std::cerr << "[Core::Genotyping] Computing backward probabilities." << std::endl;
	compute_backward_prob();
	std::cerr << "[Core::Genotyping] Computing forward probabilities and genotype likelihoods." << std::endl;
	compute_forward_prob();
}

GenotypeHMM::~GenotypeHMM()
{
	std::cerr << "[Core::Genotyping] Deleting HMM." << std::endl;
	init(backward_pass_table, 0);
	init(hmm_columns, 0);
	init(haplotype_mapper_table, 0);
	init(transition_probabilities, 0);
}

void GenotypeHMM::clear_backward_table()
{
	size_t column_count = column_iterator.get_column_count();
	init(backward_pass_table, column_count);
}

unique_ptr<vector<uint32_t> > GenotypeHMM::extract_read_ids(const vector<const Entry *>& entries) {
	unique_ptr<vector<uint32_t> > read_ids(new vector<uint32_t>());
	for (uint32_t i=0; i < entries.size(); i++) {
		read_ids->push_back(entries[i]->get_read_id());
	}
	return read_ids;
}

void GenotypeHMM::compute_index(){
	size_t column_count = column_iterator.get_column_count();
	if(column_count == 0) return;
	init(hmm_columns, column_count);
	init(haplotype_mapper_table, column_count);
	init(transition_probabilities, column_count-1);
	// do one forward pass to get the indexers (that are needed in forward and backward pass)
	column_iterator.jump_to_column(0);
	unique_ptr<vector<const Entry*> > current_input_column;
	unique_ptr<vector<const Entry*> > next_input_column;
	unique_ptr<vector<uint32_t> > current_read_ids;
	unique_ptr<vector<uint32_t> > next_read_ids;
	Column* current_column = nullptr;
	next_input_column = column_iterator.get_next();
	next_read_ids = extract_read_ids(*next_input_column);
	long double transition_constant = 0.000004L * ((long double)recombrate) * ((long double)eff_pop_size);

	for(size_t column_index=0; column_index < column_iterator.get_column_count(); ++column_index){
		
		current_input_column = std::move(next_input_column);
		current_read_ids = std::move(next_read_ids);
		if (column_iterator.has_next()) {
			next_input_column = column_iterator.get_next();
			next_read_ids = extract_read_ids(*next_input_column);
			current_column = new Column(*current_read_ids, *next_read_ids, read_set);
			hmm_columns[column_index] = current_column;
			haplotype_mapper_table[column_index] = new HaplotypeMapper(variant_info_table->at(column_index).genotype_likelihoods, variant_info_table->at(column_index).allele_references);
			transition_probabilities[column_index] = new TransitionProbabilities(calculate_transition_probabilities(variant_info_table->at(column_index).position, variant_info_table->at(column_index+1).position, transition_constant, num_haplotypes));
		} 
		else {
			assert (column_index == column_iterator.get_column_count() - 1);
			current_column = new Column(*current_read_ids, vector<uint32_t>{}, read_set); 
			hmm_columns[column_index] = current_column;
			haplotype_mapper_table[column_index] = new HaplotypeMapper(variant_info_table->at(column_index).genotype_likelihoods, variant_info_table->at(column_index).allele_references);
		}
	}
}

void GenotypeHMM::compute_backward_prob()
{
	clear_backward_table();
	uint32_t column_count = column_iterator.get_column_count();

	// set active column to the rightmost column
	column_iterator.jump_to_column(column_count-1);
	// backward pass: create sparse table
	size_t k = (size_t)sqrt(column_count);
	for(uint32_t column_index = column_count-1; column_index >= 0; --column_index){
		compute_backward_column(column_index);
		/**
		 * To conserve space, we only keep every k columns' backward values
		 * k = (size_t)sqrt(column_count)
		 * 
		 * Later when we have to calculate genotype likelihoods and need the backward values,
		 * we identify the nearest column where the backward values are available
		 * and recalculate the values.
		 * 
		 * Eg. for 10 columns, we have k = 3
		 * The following columns are stored:
		 * Idx ->    0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9
		 * Stored->  N | N | Y | N | N | Y | N | N | Y | Y
		 */
		if ((k>1) && (column_index < column_count-1) && (((column_index+1)%k) != 0)) {
			delete backward_pass_table[column_index+1];
			backward_pass_table[column_index+1] = nullptr;
		}
	}
}

void GenotypeHMM::compute_forward_prob() {

	// reset active column to the leftmost column
	column_iterator.jump_to_column(0);
	for (size_t column_index=0; column_index < column_iterator.get_column_count(); ++column_index) {
		compute_forward_column(column_index);
	}
}

/**
 * Note:
 * In this function, the following notation will be used to refer to columns.
 * 
 * Column at index column_index -> previous/prev
 * Column at index column_index - 1 -> current/curr
 * 
 * Reason: Since this is backward pass, we have already calculated the values at column_index
 * and now calculate the values at column_index-1.
 */
void GenotypeHMM::compute_backward_column(size_t column_index) {

	// IMPORTANT: The backward_pass_column_table[column_index - 1] is filled and not backward_pass_column_table[column_index].
	//            It uses backward_pass_column_table[column_index] to calculate the next column!

	// NOTE: Need column_index = 0 since we need to store the scaling parameter for the column.
	assert(column_index < column_iterator.get_column_count());
	// if current input column was not provided, create it
	unique_ptr<vector<const Entry*>> current_input_column = nullptr;
	column_iterator.jump_to_column(column_index);
	current_input_column = column_iterator.get_prev();
	
	if(column_index > 0){
		/**
		 * checks if the current column is already filled with scores.
		 * if yes, then no need to do anything.
		 */
		if (backward_pass_table[column_index-1] != nullptr) return;
	}
	
	/**
	 * Initializing objects and retrieving appropriate information
	 */
	Column* prev_indexer;
	prev_indexer = hmm_columns[column_index];
	assert(prev_indexer != nullptr);

	vector<int> prev_haplotype_to_allele = variant_info_table->at(column_index).allele_references;     // This contains the haplotype-to-allele mapping for the position column_index		
	HaplotypeMapper* prev_haplotype_mapper = haplotype_mapper_table.at(column_index);
	uint32_t num_prev_ref_states = prev_haplotype_mapper->get_num_states();
	vector<long double>* current_backward_scores = nullptr;
	uint32_t n_alleles = variant_info_table->at(column_index).get_num_alleles();
	/**
	 * TODO: get rid of useless alleles from EmissionProbabilityComputer
	 */
	EmissionProbabilityComputer emission_probability_computer = EmissionProbabilityComputer(n_alleles);
	
	/**
	 * Declaration of variables that are needed for the calculations
	 */
	long double scaling_sum = 0.0L;	// normalization factor for the backward prob.
	uint32_t read_cluster_bit_rep;	// to store the bipartition of the read clusters
	uint32_t bipartition_index;		// to store the numeric value of the bipartition without considering constraint positions
	uint32_t r_index;				// to store the index of states inside the bipartition
	uint32_t state_index;			// combining bipartition_index and haplotype_index to get the index of the particular state.
	uint32_t prev_allele0;    		// to store allele0 of previous column
	uint32_t prev_allele1;			// to store allele1 of previous column
	uint32_t prev_r_index;        	// to store the index inside bipartition of previous column
	pair<uint32_t, uint32_t> haplotypes;	// storing haplotypes
	vector<uint32_t> compatible_bipartitions; // to store bipartition indices of current column
	
	/**
	 * Helper variables
	 */
	long double beta_helper_0;		 // This helper just store the value of beta(R1,R2).
	long double beta_helper_1;       // This helper value is the beta(*,*) value.
	vector<long double> beta_helper_2(num_haplotypes);      // This helper value is the beta(R1,*) value.
	vector<long double> beta_helper_3(num_haplotypes);      // This helper value is the beta(*,R2) value.

	vector<long double>* previous_backward_scores = nullptr;
	/**
	 * Check if there are backward scores from previous column.
	 * If we are at column_count - 1, then they don't exist and we initialize it
	 */
	if(column_index < column_iterator.get_column_count()-1){
		previous_backward_scores = backward_pass_table[column_index];
	}
	else {
		backward_pass_table[column_index] = new vector<long double>(prev_indexer->get_num_bipartition()*num_prev_ref_states, 1.0L);
		previous_backward_scores = backward_pass_table[column_index];
	}
	// calculating variables required for all columns other than column 0 (initilization column)
	if (column_index > 0) {	
		/**
		 * Computing the transition probabilities of Li Stephens model
		 */
		TransitionProbabilities* transition_probability = transition_probabilities[column_index-1];
		HaplotypeMapper* curr_haplotype_mapper = haplotype_mapper_table.at(column_index-1);
		uint32_t num_curr_ref_states = curr_haplotype_mapper->get_num_states();
		Column* curr_indexer = hmm_columns.at(column_index-1);
		unique_ptr<BipartitionIterator> iterator = prev_indexer->get_iterator(read_set);
		/**
		 * determine the number of elements needed at this position
		 * number of states = number of bipartitions * number of states per bipartition
		 *                           |                                      |
		 *                           V                                      V
		 *              from the BipartitionIterator              from HaplotypeMapper
		 */
		uint32_t num_total_states = curr_indexer->get_num_bipartition() * num_curr_ref_states;
		current_backward_scores = new vector<long double>(num_total_states, 0.0L);
		vector<int> curr_haplotype_to_allele = variant_info_table->at(column_index-1).allele_references;     // This contains the haplotype-to-allele mapping for the position column_index-1

		/**
		 * This iterator iterates through all the bipartitions of previous column.
		 */
		while (iterator->has_next()){
			int bit_changed = -1;
			iterator->advance(&bit_changed);
			// Update the emission probability based on the bipartition defined by the iterator
			emission_probability_computer.update_emission_probability(bit_changed, *iterator, *current_input_column);
			/**
			 * getting the indices from the iterator.
			 * bipartition_index gives the bipartition number as determined by the Gray Code.
			 * read_cluster_bit_rep gives the info of which read clusters are in which bipartition.
			 */
			bipartition_index = iterator->get_bipartition_index();
			read_cluster_bit_rep = iterator->get_read_cluster_bit_representation();

			/**
			 * Here we use curr_indexer = hmm_columns[column_index - 1] (prev_indexer is hmm_columns[column_index]) 
			 * since we want to find the compatible bipartition is column_index - 1 (the column for which to 
			 * calculate backward probabilities).
			 * 
			 * This cant be done with prev_indexer since it has columns column_index and column_index + 1.
			 * 
			 * Hence we use curr_indexer which has columns column_index - 1 and column_index.
			 */
			curr_indexer->get_backward_compatible_bipartitions(read_cluster_bit_rep, compatible_bipartitions);
			
			/**
			 * Resetting helper variables
			 */
			beta_helper_1 = 0.0L;
			beta_helper_2.assign(num_haplotypes, 0.0L);
			beta_helper_3.assign(num_haplotypes, 0.0L);

			/**
			 * Iterating over all the states inside a bipartition given by HaplotypeMapper
			 *   at column_index.
			 * We construct the helper variables by looking at the map from
			 * index in the reduced linear-space vector to the haplotype pair.
			 */
			for (r_index = 0; r_index < prev_haplotype_mapper->get_num_states(); r_index++) {
				haplotypes = prev_haplotype_mapper->get_haplotypes_indices(r_index);
				
				// Calculating beta x emission for the helper variables.
				state_index = get_node_index(bipartition_index, r_index, num_prev_ref_states);
				prev_allele0 = prev_haplotype_to_allele[haplotypes.first];
				prev_allele1 = prev_haplotype_to_allele[haplotypes.second];
				long double b = previous_backward_scores->at(state_index) * emission_probability_computer.at(prev_allele0, prev_allele1);
				// Updating the helper variables
				beta_helper_1 += b;
				beta_helper_2[haplotypes.first] += b;
				beta_helper_3[haplotypes.second] += b;
				scaling_sum += previous_backward_scores->at(state_index);       // Updating scaling sum for later normalization
			}
			/**
			 * We have calculated the helper variables for a bipartition in the previous column (Bx).
			 * Now we find all the bipartitions in the current column (By) which are compatible with 
			 *    the bipartition of the previous column.
			 * We calculate the contribution of Bx on the different By.
			 */
			for (uint32_t compatible_bipartition_index: compatible_bipartitions) {
				/**
				 * Now we are iterating over the states in some bipartition
				 * of the CURRENT column.
				 */
				for (r_index = 0; r_index < num_curr_ref_states; r_index++) {
					/**
					 * This is the index for current column at compatible bipartition given by compatible_bipartition_index
					 * index inside the bipartition given by r_index
					 */
					state_index = get_node_index(compatible_bipartition_index, r_index, num_curr_ref_states);
					/**
					 * We find what haplotypes the r_index is pointing to.
					 * Let's say it is hap_i and hap_j.
					 * 
					 * Now we find out which alleles was present in the PREVIOUS column
					 * with these haplotpyes i and j.
					 * 
					 * We do this to calculate beta_helper_0 which is the backward score of 
					 * the state in the previous column with the same haplotypes (multiplied with
					 * the emission of that state).
					 */
					haplotypes = curr_haplotype_mapper->get_haplotypes_indices(r_index);
					if (prev_haplotype_mapper->get_state_index(haplotypes.first, haplotypes.second) == -1) {
						/**
						 * This indicates that the hap_i and hap_j pair we are looking at in the current column
						 * does not exist in the previous column.
						 * 
						 * This is most likely because those haplotypes in the previous position had a genotype
						 * which was not selected.
						 * 
						 * Since the state does not exist, the beta value is 0.
						 */
						beta_helper_0 = 0.0L;
					}
					else {
						/**
						 * The haplotype pairs have a valid entry in the previous column.
						 * 
						 * Extracting the alleles the haplotypes refer to, and the index inside
						 * the bipartitions at previous column.
						 */
						prev_allele0 = prev_haplotype_to_allele[haplotypes.first];
						prev_allele1 = prev_haplotype_to_allele[haplotypes.second];
						prev_r_index = prev_haplotype_mapper->get_state_index(haplotypes.first, haplotypes.second);
						/**
						 * the index of the state in the previous column is given by bipartition_index (from the iterator),
						 * the prev_r_index determined from the haplotypes of r_index, and the size of each bipartition in
						 * the previous column
						 */
						beta_helper_0 = 
							previous_backward_scores->at(get_node_index(bipartition_index, prev_r_index, num_prev_ref_states)) * 
							emission_probability_computer.at(prev_allele0, prev_allele1);
					}
					current_backward_scores->at(state_index) += 
						(transition_probability->q2 * beta_helper_0)
						+ (transition_probability->pq * (beta_helper_2[haplotypes.first] + beta_helper_3[haplotypes.second] - (2 * beta_helper_0)))
						+ (transition_probability->p2 * (beta_helper_1 - beta_helper_2[haplotypes.first] - beta_helper_3[haplotypes.second] + beta_helper_0));
				}
			}
		}
	}
	else {
		/**
		 * This block is executed when we are at column_index = 0
		 * So we have finished calculating all the backward values for each column
		 *  (since at column_index = 1, we did the calculation for 0)
		 * 
		 * But the values from the calculation of column_index 0 still need to be normalized.
		 * 
		 * So here we calculating the normalizing factor.
		 */
		 for (uint32_t index = 0; index < previous_backward_scores->size(); index++) scaling_sum += previous_backward_scores->at(index);
	}
	/**
	 * Now we normalize backward scores of column_index.
	 * 
	 * Here we see why the compute_backward_column() was executed for column_index = 0.
	 * Here we normalize that.
	 */
	if(previous_backward_scores != nullptr){
		std::transform((*previous_backward_scores).begin(), (*previous_backward_scores).end(), (*previous_backward_scores).begin(), [scaling_sum](long double val) { return val/scaling_sum; });
	}
	/**
	 * We also divide all the backward values calculated for column_index - 1 with the same factor.
	 * This does not normalize the values but makes it less likely that we run into underflow issues
	 *   when we use values of column_index - 1 to calculate for column_index - 2.
	 */
	if(current_backward_scores != nullptr){
		std::transform((*current_backward_scores).begin(), (*current_backward_scores).end(), (*current_backward_scores).begin(), [scaling_sum](long double val) { return val/scaling_sum; });
		backward_pass_table[column_index-1] = current_backward_scores;
	}
	scaling_parameters[column_index] = scaling_sum;
}

/**
 * Note:
 * In this function, the following notation will be used to refer to columns.
 * 
 * Column at index column_index -> current/curr
 * Column at index column_index - 1 -> previous/prev
 * 
 * Reason: Since this is forward pass, we have already calculated the values at column_index - 1
 * and now calculate the values at column_index.
 */
void GenotypeHMM::compute_forward_column(size_t column_index)
{
	assert(column_index < column_iterator.get_column_count());

	/**
	 * To conserve space, we have stored only some columns' backward values.
	 * Now that we want to re-calculate the backward values for this column if it was not stored.
	 * 
	 * How is this done?
	 *  - We find that nearest column where the values are stored.
	 *  - Use that as a starting point to re-compute the backward values until our column.
	 * 
	 * Building on the example from compute_backward_prob()
	 * for 10 columns, we have k = 3
	 * The following columns are stored:
	 * Idx ->    0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9
	 * BP->      N | N | Y | N | N | Y | N | N | Y | Y
	 * 
	 * We have to calculate GL (Genotype Likelihoods) for idx 0. But no BP (Backward Prob) available.
	 * Closest idx where values are known is idx = 2. So from idx = 2, we calculate BP for idx = 1 and idx = 0.
	 * 
	 * Idx ->    0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9
	 * BP->      Y | Y | Y | N | N | Y | N | N | Y | Y
	 * FP->      Y | N | N | N | N | N | N | N | N | N
	 * GL->      Y | N | M | N | N | N | N | N | N | N
	 * 
	 * Now we have to calculate GL for idx = 1. We have BP for idx = 1!
	 * Even though it was not originally stored, it was re-computed and kept!
	 * 
	 * Idx ->    0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9
	 * BP->      N | Y | Y | N | N | Y | N | N | Y | Y
	 * FP->      N | Y | N | N | N | N | N | N | N | N
	 * GL->      Y | Y | M | N | N | N | N | N | N | N
	 * 
	 * Note: once we have calculated GL for idx = 0, we delete the BP and FP at idx = 0.
	 * 
	 * Highlights of this strategy:
	 *  - We store a max of O(k) columns at any point of time. (Reduced space complexity)
	 * 	- Each column (on average) requires 2 Backward Pass (for initial compute and re-compute) and 1 Forward Pass.
	 *    (We increase time but the time complexity remains same)
	 */
	size_t k = (size_t)sqrt(column_iterator.get_column_count());
	vector<long double>* backward_probabilities = nullptr;
	backward_probabilities = backward_pass_table[column_index];
	// if column is not stored, recompute it
	if(backward_probabilities == nullptr) {
		// compute index of next column that has been stored
		size_t next = std::min((uint32_t) ( ((column_index + k) / k) * k ), column_iterator.get_column_count()-1);
		for(size_t i = next; i > column_index; --i){
			compute_backward_column(i);
		}
		if (backward_pass_table[column_index] ==  nullptr) {
			compute_backward_column(next);
			delete backward_pass_table[next-1];
			backward_pass_table[next-1] = nullptr;
		}
		assert (backward_pass_table[column_index] != nullptr);
		// last column just computed still needs to be scaled
		long double scaling_sum = scaling_parameters[column_index];
		std::transform((*backward_pass_table[column_index]).begin(), (*backward_pass_table[column_index]).end(), (*backward_pass_table[column_index]).begin(), [scaling_sum](long double val) { return val/scaling_sum; });
	}
	backward_probabilities = backward_pass_table[column_index];
	assert(backward_probabilities != nullptr);

	// Get the active entries at this position
	unique_ptr<vector<const Entry*>> current_input_column = nullptr;
	column_iterator.jump_to_column(column_index);
	current_input_column = column_iterator.get_next();
	
	/**
	 * Initializing objects and retrieving appropriate information
	 */
	Column* curr_indexer = hmm_columns[column_index];
	assert(curr_indexer != nullptr);
	uint32_t num_curr_bipartitions = curr_indexer->get_num_bipartition();
	Column* prev_indexer;
	uint32_t num_prev_bipartitions;
	HaplotypeMapper* curr_haplotype_mapper = haplotype_mapper_table.at(column_index);
	uint32_t num_curr_ref_states = curr_haplotype_mapper->get_num_states();
	uint32_t n_alleles = variant_info_table->at(column_index).get_num_alleles();
	/**
	 * TODO: get rid of useless alleles from EmissionProbabilityComputer
	 */
	EmissionProbabilityComputer emission_probability_computer = EmissionProbabilityComputer(n_alleles);
	HaplotypeMapper* prev_haplotype_mapper;
	uint32_t num_prev_ref_states;
		

	/**
	 * Declaration of variables
	 */
	// sum of alpha*beta, used to normalize the likelihoods
	long double normalization = 0.0L;
	/**
	 * sum of all alpha values to normalize the forward values
	 * to prevent underflow errors in next column's forward values.
	 */
	long double sum = 0.0L;
	uint32_t read_cluster_bit_rep;	// to store the bipartition of the read clusters
	uint32_t bipartition_index;		// to store the numeric value of the bipartition without considering constraint positions
	uint32_t r_index;				// to store the index of states inside the bipartition
	uint32_t state_index;			// combining bipartition_index and haplotype_index to get the index of the particular state.
	vector<uint32_t> compatible_bipartitions; // to store bipartition indices of current column
	vector<int> curr_haplotype_to_allele = variant_info_table->at(column_index).allele_references;     // This contains the haplotype-to-allele mapping for the position column_index
	vector<int> prev_haplotype_to_allele;
	pair<uint32_t, uint32_t> haplotypes;	// storing haplotypes
	TransitionProbabilities* transition_probability;
	uint32_t curr_allele0;    		// to store allele0 of current column
	uint32_t curr_allele1;			// to store allele1 of current column
	uint32_t prev_r_index;			// to store r_index of the previous column.
	
	/**
	 * Resizing helpers from current column
	 * This will be used for calculating values for column_index + 1 (the column after current column)
	 */
	curr_alpha_helper_1.assign(num_curr_bipartitions, 0.0L);
	curr_alpha_helper_2.resize(num_curr_bipartitions);
	for (auto& row : curr_alpha_helper_2) {
		row.assign(num_haplotypes, 0.0L);
	}
	curr_alpha_helper_3.resize(num_curr_bipartitions);
	for (auto& row : curr_alpha_helper_3) {
		row.assign(num_haplotypes, 0.0L);
	}
	
	/**
	 * Helper variables to get the values from previous column
	 */
	long double ah_0;		// takes the value of alpha(R1, R2)
	long double ah_1;       // Takes the value of alpha_helper_1->at(bipartition_index)
	vector<long double> ah_2;      // Takes the value of alpha_helper_2->at(bipartition_index)
	vector<long double> ah_3;      // Takes the value of alpha_helper_3->at(bipartition_index)

	/**
	 * initializing the vector to store forward probabilities of current column
	 */
	uint32_t num_total_states = num_curr_bipartitions * num_curr_ref_states;
	current_forward_probabilities.assign(num_total_states, 0.0L);
	

	// calculating variables required for all columns other than column 0 (initilization column)
	if (column_index > 0) {
		/**
		 * Computing the transition probabilities of Li Stephens model
		 */
		transition_probability = transition_probabilities[column_index-1];
		prev_haplotype_to_allele = variant_info_table->at(column_index-1).allele_references;     // This contains the haplotype-to-allele mapping for the position column_index-1
		prev_indexer = hmm_columns[column_index-1];
		num_prev_bipartitions = prev_indexer->get_num_bipartition();
		prev_haplotype_mapper = haplotype_mapper_table.at(column_index-1);
		num_prev_ref_states = prev_haplotype_mapper->get_num_states();
	}
	// iterate over all bipartitions
	unique_ptr<BipartitionIterator> iterator = curr_indexer->get_iterator(read_set);
	while (iterator->has_next()) {
		int bit_changed = -1;
		iterator->advance(&bit_changed);
		// Update the emission probability based on the bipartition defined by the iterator
		emission_probability_computer.update_emission_probability(bit_changed, *iterator, *current_input_column);
		/**
		 * getting the indices from the iterator.
		 * bipartition_index gives the bipartition number as determined by the Gray Code.
		 * read_cluster_bit_rep gives the info of which read clusters are in which bipartition.
		 */
		bipartition_index = iterator->get_bipartition_index();
		read_cluster_bit_rep = iterator->get_read_cluster_bit_representation();
		if (column_index == 0) {
			/**
			 * Calculating the forward probabilities for the first column.
			 * Since no prior information is available, this column gets the values of the emissions.
			 */
			for (r_index = 0; r_index < num_curr_ref_states; r_index++) {
				haplotypes = curr_haplotype_mapper->get_haplotypes_indices(r_index);
				state_index = get_node_index(bipartition_index, r_index, num_curr_ref_states);
				curr_allele0 = curr_haplotype_to_allele[haplotypes.first];
				curr_allele1 = curr_haplotype_to_allele[haplotypes.second];
				// storing the emissions in the forward probability vector.
				current_forward_probabilities[state_index] = emission_probability_computer.at(curr_allele0, curr_allele1); 
			}
		}
		else {
			/**
			 * Here we use prev_indexer = hmm_columns[column_index - 1] (curr_indexer is hmm_columns[column_index]) 
			 * since we want to find the compatible bipartition is column_index - 1
			 * 
			 * This cant be done with curr_indexer since it has columns column_index and column_index + 1.
			 * 
			 * Hence we use prev_indexer which has columns column_index - 1 and column_index.
			 */
			/**
			 * Iterating through compatible bipartitions of previous column
			 */
			prev_indexer->get_backward_compatible_bipartitions(read_cluster_bit_rep, compatible_bipartitions);
			for (uint32_t compatible_bipartition_index : compatible_bipartitions) {
				assert(compatible_bipartition_index < num_prev_bipartitions);
				/**
				 * Gathering the helpers for this bipartition
				 */
				ah_1 = alpha_helper_1[compatible_bipartition_index];
				ah_2 = alpha_helper_2[compatible_bipartition_index];
				ah_3 = alpha_helper_3[compatible_bipartition_index];
				/**
				 * Iterating through all the states inside bipartition_index of current column.
				 * In this loop, we will calculate the contributions of the states in compatible_bipartition_index (of previous column)
				 * to the states in bipartition_index (of current column).
				 * 
				 * For this calculation, we will use the helper variable of compatible_bipartition_index.
				 */
				for (r_index = 0; r_index < num_curr_ref_states; r_index++) {
					/**
					 * We find what haplotypes the r_index is pointing to.
					 * Let's say it is hap_i and hap_j.
					 * 
					 * Now we find out which alleles was present in the PREVIOUS column
					 * with these haplotpyes i and j.
					 * 
					 * We do this to calculate ah_0 which is the forward score of the state in
					 * previous column which has hap_i and hap_j and compatible_bipartition_index as bipartition.
					 */
					haplotypes = curr_haplotype_mapper->get_haplotypes_indices(r_index);
					curr_allele0 = curr_haplotype_to_allele.at(haplotypes.first);
					curr_allele1 = curr_haplotype_to_allele.at(haplotypes.second);
					/**
					 * This is the index for current column at compatible bipartition given by compatible_bipartition_index
					 * index inside the bipartition given by r_index
					 */
					state_index = get_node_index(bipartition_index, r_index, num_curr_ref_states);
					if (prev_haplotype_mapper->get_state_index(haplotypes.first, haplotypes.second) == -1) {
						/**
						 * This indicates that the hap_i and hap_j pair we are looking at in the current column
						 * does not exist in the previous column.
						 * 
						 * This is most likely because those haplotypes in the previous position had a genotype
						 * which was not selected.
						 * 
						 * Since the state does not exist, the ah_0 is 0.
						 */
						ah_0 = 0.0L;
					}
					else {
						/**
						 * The haplotype pairs have a valid entry in the previous column.
						 * 
						 * ah_0 will get the forward score value of state given by compatible_bipartition_index
						 * and hap_i and hap_j.
						 */
						prev_r_index = prev_haplotype_mapper->get_state_index(haplotypes.first, haplotypes.second);
						ah_0 = previous_forward_probabilities[get_node_index(compatible_bipartition_index, prev_r_index, num_prev_ref_states)];
						
					}
					// Updating the forward value
					current_forward_probabilities[state_index] += 
						emission_probability_computer.at(curr_allele0, curr_allele1)
							* ((transition_probability->q2 * ah_0 )
							+ (transition_probability->pq * (ah_2[haplotypes.first] + ah_3[haplotypes.second] - (2 * ah_0)))
							+ (transition_probability->p2 * (ah_1 - ah_2[haplotypes.first] - ah_3[haplotypes.second] + ah_0)));
					/**
					 * Updating the alpha helpers of current column 
					 */	
					curr_alpha_helper_1[bipartition_index] += current_forward_probabilities[state_index];
					curr_alpha_helper_2[bipartition_index][haplotypes.first] += current_forward_probabilities[state_index];
					curr_alpha_helper_3[bipartition_index][haplotypes.second] += current_forward_probabilities[state_index];
				}
			}
		}
	}
	
	/**
	 * Calculating the Genotype Likelihoods of current column
	 */
	long double forward_backward = 0.0L;
	uint32_t cannonical_genotype_index;
	vector<uint32_t> sorted_alleles;
	sorted_alleles.reserve(2);
	variant_information_t variant_info = variant_info_table->at(column_index);
	assert (current_forward_probabilities.size() == backward_probabilities->size());
	variant_info.genotype_likelihoods.reset(); // reseting the likelihood vector since it still has values from last genotyping round.
	for (bipartition_index = 0; bipartition_index < num_curr_bipartitions; bipartition_index++) {
		for (r_index = 0; r_index < num_curr_ref_states; r_index++) {
			state_index = get_node_index(bipartition_index, r_index, num_curr_ref_states);
			haplotypes = curr_haplotype_mapper->get_haplotypes_indices(r_index);
			curr_allele0 = curr_haplotype_to_allele[haplotypes.first];
			curr_allele1 = curr_haplotype_to_allele[haplotypes.second];
			sorted_alleles.clear();
			if (curr_allele0 < curr_allele1) {
				sorted_alleles.push_back(curr_allele0);
				sorted_alleles.push_back(curr_allele1);
			}
			else {
				sorted_alleles.push_back(curr_allele1);
				sorted_alleles.push_back(curr_allele0);
			}
			cannonical_genotype_index = convert_alleles_to_index(sorted_alleles);
			current_forward_probabilities[state_index] /= scaling_parameters[column_index];
			sum += current_forward_probabilities[state_index];
			forward_backward = current_forward_probabilities[state_index] * backward_probabilities->at(state_index);
			normalization += forward_backward;

			variant_info.genotype_likelihoods.increment_by_index(cannonical_genotype_index, forward_backward);
		}
	}
	
	// normalzing the forward probabilities of current column
	std::transform(current_forward_probabilities.begin(), current_forward_probabilities.end(), current_forward_probabilities.begin(), [sum](long double val) { return val/sum; });
	// normalize the helper variables
	std::transform(curr_alpha_helper_1.begin(), curr_alpha_helper_1.end(), curr_alpha_helper_1.begin(), [sum](long double val) { return val/sum; });
	for (uint32_t i = 0; i < num_curr_bipartitions; i++) {
		std::transform(curr_alpha_helper_2[i].begin(), curr_alpha_helper_2[i].end(), curr_alpha_helper_2[i].begin(), [sum](long double val) { return val/sum; });
		std::transform(curr_alpha_helper_3[i].begin(), curr_alpha_helper_3[i].end(), curr_alpha_helper_3[i].begin(), [sum](long double val) { return val/sum; });
	}
	// normalize the likelihoods
	variant_info.genotype_likelihoods.divide_likelihoods_by(normalization);

	// update the variant info tables active alleles based on the calculated likelihoods
	std::vector<uint32_t> selected_genotype_indices = variant_info.genotype_likelihoods.select_genotypes();
	variant_info.update_active_alleles(ploidy, selected_genotype_indices);
	if (variant_info.phasable) {
		read_set->setEntryAlleles(variant_info.position, variant_info.active_alleles);
	}

	/**
	 * Replace the forward values from previous column to current column.
	 * 
	 * Using a swap function to avoid creation and destruction of memory space.
	 */
	swap(previous_forward_probabilities, current_forward_probabilities);
	
	/**
	 * Replacing the alpha helpers from previous column with
	 * the helpers calculated with current column.
	 * 
	 * Using a swap function to avoid creation and destruction of memory space.
	 */
	swap(alpha_helper_1, curr_alpha_helper_1);
	swap(alpha_helper_2, curr_alpha_helper_2);
	swap(alpha_helper_3, curr_alpha_helper_3);
}

vector<long double> GenotypeHMM::get_genotype_likelihoods(uint32_t index) {
	assert(index < column_iterator.get_column_count());
	return variant_info_table->at(index).genotype_likelihoods.as_vector();
}


uint32_t GenotypeHMM::get_node_index(uint32_t b_index, uint32_t r_index, uint32_t num_states) {
	return (b_index * num_states) + r_index;
}