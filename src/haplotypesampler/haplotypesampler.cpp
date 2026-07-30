/**
 * Original implementation at https://github.com/eblerjana/pangenie/blob/master/src/haplotypesampler.cpp (Commit f682fb6)
 */

#include "haplotypesampler.h"
#include "../columniterator.h"
#include "samplingtransitions.h"
#include <algorithm>
#include <cassert>
#include <map>
#include <math.h>
#include <fstream>

using namespace std;

void print_dpcolumn(DPColumn* column) {
	cout << "Print column:" << endl;
	for (uint32_t i = 0; i < column->column.size(); ++i) {
		cout << column->column.at(i) << endl;
	}
	cout << "--------" << endl;
}

HaplotypeSampler::HaplotypeSampler(
	ReadSet* read_set,
	std::vector<variant_information_t>* variant_info_table,	
	const uint32_t size,
	const float recombrate, 
	const float effective_N, 
	vector<uint32_t>* best_scores, 
	const bool remove_reference, 
	const uint32_t allele_penalty)
	: read_set(read_set),
	variant_info_table(variant_info_table),
	recombrate(recombrate),
	effective_N(effective_N),
	allele_penalty(allele_penalty) {
	
	if (size < 1) return;

	compute_emission_probabilities();
	num_haplotypes = variant_info_table->at(0).allele_references.size();
	
	// generate size Viterbi paths
	for (uint32_t i = 0; i < size; ++i) {
		compute_viterbi_path(best_scores);
	}

	if (!remove_reference) this->sampled_paths.sampled_paths.push_back(vector<uint32_t>(variant_info_table->size(), 0));	
	// print the sampled paths if requested
	/*
	if (path_output != "") {
		ofstream path_outfile;
		path_outfile.open(path_output);

		// print the header line
		path_outfile << "#chromosome\tposition";
		for (uint32_t path_id = 0; path_id < this->sampled_paths.sampled_paths.size(); ++path_id) {
			path_outfile << "\tHaplotypeID_path" << path_id << "\tRecombination_path" << path_id;
		}
		path_outfile << endl;

		// print stats for each position
		for (uint32_t column_index = 0; column_index < unique_kmers->size(); ++column_index) {
			path_outfile << chromosome << "\t" << this->unique_kmers->at(column_index)->get_variant_position();
			for (uint32_t path_id = 0; path_id < this->sampled_paths.sampled_paths.size(); ++path_id) {
				path_outfile << "\t" << this->sampled_paths.sampled_paths.at(path_id).at(column_index);
				path_outfile << "\t" << this->sampled_paths.recombination(column_index, path_id);
			}
			path_outfile << endl;
		}
	}
	*/

	// clean up
	init(this->viterbi_columns, 0);
	init(this->viterbi_backtrace_columns, 0);
	
}

void HaplotypeSampler::compute_emission_probabilities() {

	this->emission_costs.reserve(variant_info_table->size());
	
	ColumnIterator iterator(*read_set, variant_info_table);
	unique_ptr<vector<const Entry*> > current_input_column;
	iterator.jump_to_column(0);
	uint32_t column_count = iterator.get_column_count();
	assert (column_count == variant_info_table->size());
	
	for (uint32_t column_index = 0; column_index < column_count; ++column_index) {
		assert(iterator.has_next());
		current_input_column = iterator.get_next();
		this->emission_costs.push_back(SamplingEmissions(*current_input_column, variant_info_table->at(column_index).get_num_alleles()));	
	}
	assert(!iterator.has_next());
}


void HaplotypeSampler::get_column_minima(std::vector<uint32_t>& column, std::vector<bool>& mask, uint32_t& first_id, uint32_t& second_id, uint32_t& first_val, uint32_t& second_val) const
{
 
	assert (column.size() > 1);
	assert (column.size() == mask.size());
 
    first_val = std::numeric_limits<uint32_t>::max();
	second_val = std::numeric_limits<uint32_t>::max();

	first_id = std::numeric_limits<uint32_t>::max();
	second_id = std::numeric_limits<uint32_t>::max();

    for (uint32_t i = 0; i < column.size(); i++) {
		// if masked with false, ignore.
		if (!mask[i]) continue;

        /* If current element is smaller than first
        then update both first and second */
        if (column[i] < first_val) {
            second_val = first_val;
			second_id = first_id;
			first_val = column[i];
            first_id = i;
        } else if ( (column[i] < second_val) && (i != first_id)) {
            second_val = column[i];
			second_id = i;
    	}
	}
}


void HaplotypeSampler::compute_viterbi_path(vector<uint32_t>* best_scores) {
	uint32_t column_count = this->variant_info_table->size();
	init(this->viterbi_columns, column_count);
	init(this->viterbi_backtrace_columns, column_count);

	// perform Viterbi algorithm
	uint32_t k = (uint32_t) sqrt(column_count);
	for (uint32_t column_index = 0; column_index < column_count; ++column_index) {
		compute_viterbi_column(column_index);
		// store sparse table. Check if previous column needs to be deleted.
		if ((k > 1) && (column_index > 0) && (((column_index - 1)%k != 0)) ) {
			delete this->viterbi_columns[column_index-1];
			this->viterbi_columns[column_index-1] = nullptr;
			delete this->viterbi_backtrace_columns[column_index-1];
			this->viterbi_backtrace_columns[column_index-1] = nullptr;
		}
	}

	// find the best value in the last column
	DPColumn* last_column = this->viterbi_columns.at(column_count-1);
	assert (last_column != nullptr);
	uint32_t best_index = 0;
	uint32_t best_value = last_column->column.at(0);
	for (uint32_t i = 1; i < last_column->column.size(); ++i) {
		uint32_t entry = last_column->column.at(i);
		if (entry < best_value) {
			best_value = entry;
			best_index = i;
		}
	}

	// keep track of best DP score
	if (best_scores != nullptr) {
		best_scores->push_back(best_value);
	}


	// backtracking
	vector<uint32_t> path(column_count);
	uint32_t column_index = column_count - 1;
	while(true) {
		// columns might have to be re-computed
		if (this->viterbi_backtrace_columns[column_index] == nullptr) {
			uint32_t j = column_index / k*k;
			assert (this->viterbi_columns[j] != nullptr);
			for (j = j+1; j <= column_index; ++j) {
				compute_viterbi_column(j);
			}
		}
		// store the best path
		path[column_index] = best_index;

		// penalize allele covered by selected path
		unsigned short best_allele = this->variant_info_table->at(column_index).allele_references[best_index];
		emission_costs.at(column_index).penalize(best_allele, allele_penalty);

		if (column_index == 0) break;

		// update the best index
		best_index = this->viterbi_backtrace_columns.at(column_index)->at(best_index);

		// current column is no longer needed. Delete it.
		delete this->viterbi_columns[column_index];
		this->viterbi_columns[column_index] = nullptr;
		delete this->viterbi_backtrace_columns[column_index];
		this->viterbi_backtrace_columns[column_index] = nullptr;
		column_index -= 1;
	}
	this->sampled_paths.sampled_paths.push_back(path);
}

void HaplotypeSampler::compute_viterbi_column(uint32_t column_index) {
	assert (column_index < this->variant_info_table->size());

	// check whether column was computed already
	if (this->viterbi_columns[column_index] != nullptr) return;

	// get previous column
	DPColumn* previous_column = nullptr;
	vector<bool> prev_mask;

	if (column_index > 0) {
		previous_column = this->viterbi_columns[column_index - 1];

		// bitvector marking which path ids shall be ignored (removed in previous DP iterations)
		prev_mask = this->sampled_paths.mask_indexes(column_index-1, num_haplotypes-1);
	}

	DPColumn* current_column = new DPColumn();
	current_column->column = vector<uint32_t>(num_haplotypes);

	// backtrace column
	vector<uint32_t>* backtrace_column = new vector<uint32_t>(num_haplotypes, numeric_limits<uint32_t>::max());

	// precompute minima for each index in current column. helper[i] contains the value of
	// the minimum value of all positions except i in previous columns.
	vector<uint32_t> helper_val(num_haplotypes);
	vector<uint32_t> helper_id(num_haplotypes);

	// currently masked indexes (removed in previous DP iterations)
	vector<bool> cur_mask = this->sampled_paths.mask_indexes(column_index, num_haplotypes-1);
 
	SamplingTransitions* transition_cost_computer = nullptr;

	if (column_index > 0) {
		// set up SamplingTransitions
		uint32_t from_variant = this->variant_info_table->at(column_index-1).position;
		uint32_t to_variant = this->variant_info_table->at(column_index).position;
		transition_cost_computer = new SamplingTransitions(from_variant, to_variant, this->recombrate, num_haplotypes, this->effective_N);

		// compute smallest and second smallest element in previous column
		uint32_t first_id, second_id;
		uint32_t first_val, second_val;
		this->get_column_minima(previous_column->column, prev_mask, first_id, second_id, first_val, second_val);
		// fill helper vector
		for (uint32_t i = 0; i < num_haplotypes; ++i) {
			if (cur_mask[i]) {
				if (i == first_id) {
					// need to consider second smallest value from previous column
					helper_val[i] = second_val;
					helper_id[i] = second_id;
				} else {
					// need to consider minimum from previous column
					helper_val[i] = first_val;
					helper_id[i] = first_id;
				}
			} else {
				helper_val[i] = numeric_limits<uint32_t>::max();
				helper_id[i] = numeric_limits<uint32_t>::max();
			}
		}
	}

	// TODO: check for overflows (see WH code)!!

	// fill DP column based on precomputed minima
	for (uint32_t i = 0; i < num_haplotypes; ++i) {
		// if current index is masked, skip.
		if (!cur_mask[i]) {
			current_column->column[i] = numeric_limits<uint32_t>::max();
			continue;
		}
		uint32_t previous_cell = 0;
		if (column_index > 0) {
			// check of previous value exists for same path (might be masked)
			// keep track of where the minimum came from and store in backtrace table
			previous_cell = helper_val[i] + transition_cost_computer->compute_transition_cost(true);

			// check if there was an overflow
			if (previous_cell < helper_val[i]) previous_cell = numeric_limits<uint32_t>::max();

			backtrace_column->operator[](i) = helper_id[i];

			if (prev_mask[i]) {
				uint32_t same = previous_column->column.at(i) + transition_cost_computer->compute_transition_cost(false);

				// check if there was an overflow
				if (same < previous_column->column.at(i)) same = numeric_limits<uint32_t>::max();

				if (same < previous_cell) {
					previous_cell = same;
					backtrace_column->operator[](i) = i;
				}
			}
		}
		// add Emission costs
		unsigned short allele = (unsigned short)this->variant_info_table->at(column_index).allele_references[i];
		current_column->column[i] = previous_cell + emission_costs.at(column_index).get_emission_cost(allele);
		
		// check if there was an overflow
		if (current_column->column[i] < previous_cell) current_column->column[i] = numeric_limits<uint32_t>::max();
	}

	// store the column and clean up
	this->viterbi_columns.at(column_index) = current_column;
	this->viterbi_backtrace_columns.at(column_index) = backtrace_column;


	if (transition_cost_computer != nullptr) {
		delete transition_cost_computer;
	}
}


std::vector<variant_information_t> HaplotypeSampler::get_updated_variant_table(uint32_t ploidy) {
	
	std::vector<variant_information_t> updated_table(variant_info_table->size());
	variant_information_t variant;
	uint32_t position;
	uint32_t n_alleles;
	std::vector<int> allele_references(sampled_paths.sampled_paths.size());
	bool sv_flag;
	uint32_t hap_index;
	int hap_allele;
	std::vector<bool> allele_presence_vector;

	for (uint32_t i = 0; i < variant_info_table->size(); i++) {
		variant = variant_info_table->at(i);
		position = variant.position;
		n_alleles = variant.get_num_alleles();
		allele_presence_vector = std::vector<bool>(n_alleles, false);
		sv_flag = variant.is_sv;
		for (uint32_t j = 0; j < sampled_paths.sampled_paths.size(); j++) {
			hap_index = sampled_paths.sampled_paths[j][i];
			hap_allele = variant.allele_references[hap_index];
			assert(hap_allele != -1);	// should not be selecting unknown haplotype
			assert(hap_allele < n_alleles);	// should be less than the number of alleles defined at the position
			allele_references[j] = hap_allele;
			allele_presence_vector[hap_allele] = true;
		}
		updated_table[i] = variant_information_t(position, ploidy, n_alleles, allele_references, sv_flag);
		updated_table[i].active_alleles = allele_presence_vector;
		if (updated_table[i].count_active_alleles() <= 2 && !sv_flag) {
			updated_table[i].phasable = true;
		} else {
			updated_table[i].phasable = false;
		}
	}
	
	return updated_table;
}