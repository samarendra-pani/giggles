#include <stdexcept>
#include <cassert>
#include <limits>
#include <fstream>
#include <array>
#include <algorithm>
#include <cmath>
#include <vector>

#include "genotypehmm.h"
#include "columnindexingiterator.h"
#include "transitionprobabilitycomputer.h"
#include "emissionprobabilitycomputer.h"
#include "binomial.h"
#include "matrixmultiplication.h"

using namespace std;

GenotypeHMM::GenotypeHMM(ReadSet* read_set, const vector<float>& recombcost, const Pedigree* pedigree, const unsigned int& n_references, const vector<unsigned int>* positions, const vector<unsigned int>* n_allele_positions,  const vector<vector<int> >* allele_references)
    :read_set(read_set),
     recombcost(recombcost),
     pedigree(pedigree),
     input_column_iterator(*read_set, positions),
     backward_input_column_iterator(*read_set, positions),
     transition_probability_table(input_column_iterator.get_column_count() - 1,nullptr),
     scaling_parameters(input_column_iterator.get_column_count(),-1.0L),
     variant_positions(positions),
     variant_n_allele_positions(n_allele_positions),
     n_references(n_references),
     allele_references(allele_references)
{
    genotype_likelihood_table = Vector2D<genotype_likelihood_t>(pedigree->size(),input_column_iterator.get_column_count());
    assert (pedigree->size() == 1);
    for (size_t i = 0; i < pedigree->size(); i ++) {
        for (size_t j = 0; j < input_column_iterator.get_column_count(); j++) {
            genotype_likelihood_table.set(i, j, genotype_likelihood_t(binomial_coefficient(n_allele_positions->at(j)+1, n_allele_positions->at(j)-1)));
        }
    }
    read_set->reassignReadIds();
    assert(input_column_iterator.get_column_count() == backward_input_column_iterator.get_column_count());

    // translate all individual ids to individual indices
    for(size_t i = 0; i<read_set->size(); ++i)
    {
        read_sources.push_back(pedigree->id_to_index(read_set->get(i)->getSampleID()));
    }
    //compute forward and backward probabilities
    compute_index();
    compute_backward_prob();
    compute_forward_prob();
}

GenotypeHMM::~GenotypeHMM()
{
    init(forward_pass_column_table,0);
    init(backward_pass_column_table, 0);
    init(hmm_columns,0);
    init(transition_probability_table,0);
}

void GenotypeHMM::clear_forward_table()
{
    size_t column_count = input_column_iterator.get_column_count();
    init(forward_pass_column_table, 1);
}

void GenotypeHMM::clear_backward_table()
{
    size_t column_count = input_column_iterator.get_column_count();
    init(backward_pass_column_table, column_count);
}

unique_ptr<vector<unsigned int> > GenotypeHMM::extract_read_ids(const vector<const Entry *>& entries) {
    unique_ptr<vector<unsigned int> > read_ids(new vector<unsigned int>());
    for (int i=0; i < entries.size(); i++) {
        read_ids->push_back(entries[i]->get_read_id());
    }
    return read_ids;
}

void GenotypeHMM::compute_index(){
    size_t column_count = input_column_iterator.get_column_count();
    if(column_count == 0) return;
    init(hmm_columns, column_count);
    // do one forward pass to get the indexers (that are needed in forward and backward pass)
    input_column_iterator.jump_to_column(0);
    unique_ptr<vector<const Entry*> > current_input_column;
    unique_ptr<vector<const Entry*> > next_input_column;
    unique_ptr<vector<unsigned int> > current_read_ids;
    unique_ptr<vector<unsigned int> > next_read_ids;
    Column* current_column = nullptr;
    next_input_column = input_column_iterator.get_next();
    next_read_ids = extract_read_ids(*next_input_column);
    
    for(size_t column_index=0; column_index < input_column_iterator.get_column_count(); ++column_index){
        
        current_input_column = std::move(next_input_column);
        current_read_ids = std::move(next_read_ids);
        if (input_column_iterator.has_next()) {
            next_input_column = input_column_iterator.get_next();
            next_read_ids = extract_read_ids(*next_input_column);
            current_column = new Column(column_index, &n_references, *current_read_ids, *next_read_ids, read_set);
            hmm_columns[column_index] = current_column;
        } 
        else {
            assert (column_index == input_column_iterator.get_column_count() - 1);
            current_column = new Column(column_index, &n_references, *current_read_ids, vector<unsigned int>{}, read_set); 
            hmm_columns[column_index] = current_column;
        }
    }
}

void GenotypeHMM::compute_backward_prob()
{
    clear_backward_table();
    unsigned int column_count = backward_input_column_iterator.get_column_count();

    // if no reads are in the read set, nothing to do
    if(backward_input_column_iterator.get_column_count() == 0){
        return;
    }
    // do backward pass, start at rightmost column
    backward_input_column_iterator.jump_to_column(column_count-1);
    // get the next column (which is left of current one)
    unique_ptr<vector<const Entry*> > current_input_column;
    unique_ptr<vector<unsigned int> > current_read_ids;
    unique_ptr<vector<const Entry*> > next_input_column = backward_input_column_iterator.get_next();
    unique_ptr<vector<unsigned int> > next_read_ids = extract_read_ids(*next_input_column);
    // backward pass: create sparse table
    size_t k = (size_t)sqrt(column_count);
    for(int column_index = column_count-1; column_index >= 0; --column_index){
        // make former next column the current one
        current_input_column = std::move(next_input_column);
        current_read_ids = std::move(next_read_ids);
        // peek ahead and get the next column
        if (backward_input_column_iterator.has_next()){
            next_input_column = backward_input_column_iterator.get_next();
            next_read_ids = extract_read_ids(*next_input_column);
        } else {
            assert(next_input_column.get() == 0);
            assert(next_read_ids.get() == 0);
        }
        // compute the backward probabilities
        if (column_index > 0) {
            transition_probability_table[column_index - 1] = new TransitionProbabilityComputer(recombcost[column_index-1], allele_references->at(column_index));
        }
        compute_backward_column(column_index,std::move(current_input_column));
        if (column_index > 0) {
            delete transition_probability_table[column_index-1];
            transition_probability_table[column_index-1] = nullptr;
        }
        // check whether to delete the previous column
        if ((k>1) && (column_index < column_count-1) && (((column_index+1)%k) != 0)) {
            delete backward_pass_column_table[column_index+1];
            backward_pass_column_table[column_index+1] = nullptr;
        }
    }
}

void GenotypeHMM::compute_forward_prob()
{
    clear_forward_table();

    // if no reads are in read set, nothing to compute
    if (input_column_iterator.get_column_count() == 0) {
        return;
    }

    // start at leftmost column (= 0th column)
    input_column_iterator.jump_to_column(0);
    // store current and next column
    unique_ptr<vector<const Entry *> > current_input_column;
    unique_ptr<vector<const Entry *> > next_input_column;
    // get the next column ahead of time
    next_input_column = input_column_iterator.get_next();
    unique_ptr<vector<unsigned int> > next_read_ids = extract_read_ids(*next_input_column);

    // forward pass: create a sparse table, storing values at every sqrt(#columns)-th position
    for (size_t column_index=0; column_index < input_column_iterator.get_column_count(); ++column_index) {
        // make former next column the current one
        current_input_column = std::move(next_input_column);
        unique_ptr<vector<unsigned int> > current_read_ids = std::move(next_read_ids);
        // peek ahead and get the next column
        if (input_column_iterator.has_next()) {
            next_input_column = input_column_iterator.get_next();
            next_read_ids = extract_read_ids(*next_input_column);
        } else {
            assert(next_input_column.get() == 0);
            assert(next_read_ids.get() == 0);
        }
        // compute forward probabilities for the current column
        if (column_index > 0 && transition_probability_table[column_index - 1] == nullptr) {
            transition_probability_table[column_index - 1] = new TransitionProbabilityComputer(recombcost[column_index-1], allele_references->at(column_index));
        }
        compute_forward_column(column_index,std::move(current_input_column));
        if (column_index > 0) {
            delete transition_probability_table[column_index-1];
            transition_probability_table[column_index-1] = nullptr;
        }
    }
}

void GenotypeHMM::compute_backward_column(size_t column_index, unique_ptr<vector<const Entry*>> current_input_column) {

    // IMPORTANT: The backward_pass_column_table[column_index - 1] is filled and not backward_pass_column_table[column_index].
    //            It uses backward_pass_column_table[column_index] to calculate the next column!

    // NOTE: Need column_index = 0 since we need to store the scaling parameter for the column.
    assert(column_index < backward_input_column_iterator.get_column_count());
    Column* current_indexer;
    TransitionProbabilityComputer* current_transition_table;
    // check if column already exists
    if(column_index > 0){
        if (backward_pass_column_table[column_index-1] != nullptr) return;
        current_indexer = hmm_columns[column_index];
        current_transition_table = transition_probability_table[column_index-1];
        assert(current_indexer != nullptr);
        assert (current_transition_table != nullptr);
    }
    // if current input column was not provided, create it
    if(current_input_column.get() == nullptr) {
        backward_input_column_iterator.jump_to_column(column_index);
        current_input_column = backward_input_column_iterator.get_next();
    }
    vector<long double>* previous_projection_column = nullptr;

    // check if there is a projection column
    if(column_index < backward_input_column_iterator.get_column_count()-1){
        previous_projection_column = backward_pass_column_table[column_index];
    }
    else {
        backward_pass_column_table[column_index] = new vector<long double>(hmm_columns[column_index]->get_column_size(), 1.0L);
        previous_projection_column = backward_pass_column_table[column_index];
    }
    // initialize the new projection column (= current index-1)
    vector<long double>* current_projection_column = nullptr;
    if(column_index > 0){
        current_projection_column = new vector<long double>(hmm_columns[column_index-1]->get_column_size(), 0.0L);
    }
    int n_alleles = variant_n_allele_positions->at(column_index);
    Vector2D<long double> emission_probability_computer = Vector2D<long double>(n_alleles, n_alleles);
    // for scaled version of forward backward alg, keep track of the sum of backward
    long double scaling_sum = 0.0L;
    vector<unsigned int> compatible_bipartitions;
    unsigned int b_index;
    unsigned int r_index;
    unsigned int index;
    // iterate over all bipartitions of the column on the right. So we are calculating the values in column current_index - 1
    if (column_index > 0) {
        unique_ptr<ColumnIndexingIterator> iterator = current_indexer->get_iterator(read_set);
        while (iterator->has_next()){
            int bit_changed = -1;
            iterator->advance(&bit_changed);
            // Update the emission probability based on the bipartition defined by the iterator
            update_emission_probability(&emission_probability_computer, bit_changed, *iterator, *current_input_column);
            // Determine the indices that are compatible with with the bipartition at position column_index
            b_index = iterator->get_b_index();
            // Here we use hmm_columns[column_index - 1] (current_indexer is hmm_columns[column_index]) since we want to find the compatible bipartition is column_index - 1 (the column for which to calculate backward probabilities).
            // This cant be done with current_indexer since it has columns column_index and column_index + 1.
            // Hence we use hmm_columns[column_index-1] which has columns column_index - 1 and column_index.
            compatible_bipartitions = hmm_columns[column_index-1]->get_backward_compatible_bipartitions(b_index, read_set);
            long double beta_helper_1 = 0.0L;       // This helper value is the beta(*,*) value.
            vector<long double> beta_helper_2(n_references, 0.0L);      // This helper value is the beta(R1,*) value.
            vector<long double> beta_helper_3(n_references, 0.0L);      // This helper value is the beta(*,R2) value.
            
            vector<int> haplotype_to_allele_prev = allele_references->at(column_index);     // This contains the haplotype-to-allele mapping for the position column_index
            vector<int> haplotype_to_allele_curr = allele_references->at(column_index-1);     // This contains the haplotype-to-allele mapping for the position column_index-1
            for (int r_index = 0; r_index < pow(n_references,2); r_index++) {
                // Extract the ref haplotype paths and alleles from r_index
                int r = r_index;    // Making a copy of r_index to extract R1, R2.
                vector<unsigned int> ref;
                vector<int> allele;
                ref.resize(2);
                allele.resize(2);
                for (int i = 0; i < 2; i++) {
                    ref[1-i] = r%n_references;
                    allele[1-i] = haplotype_to_allele_prev[ref[1-i]];
                    r = r / n_references;
                }
                if ((allele[0] == -1) || (allele[1] == -1)) {continue;}     // Skipping nodes where the allele is not defined.
                // Calculating beta x emission for the helper variables.
                index = current_indexer->get_index(b_index, r_index);
                long double b = previous_projection_column->at(index) * emission_probability_computer.at(allele[0],allele[1]);
                // Updating the helper variables
                beta_helper_1 += b;
                beta_helper_2[ref[0]] += b;
                beta_helper_3[ref[1]] += b;
                scaling_sum += previous_projection_column->at(index);       // Updating scaling sum for later normalization
            }
            // Looping through the compatible bipartitions in column_index - 1
            for (int i = 0; i < compatible_bipartitions.size(); i++) {
                index = compatible_bipartitions[i] * pow(n_references,2);   // This is the index for column in column_index - 1. So this index corresponds to the (0,0) ref haplotype of the bipartition. So to access the rest of the indices, we need to sum index and r_index.
                // Remember that b_index is the index for column_index
                for (int r_index = 0; r_index < pow(n_references,2); r_index++) {
                    // Extract the ref haplotype paths and alleles from r_index
                    int r = r_index;    // Making a copy of r_index to extract R1, R2.
                    vector<unsigned int> ref;
                    vector<int> allele_curr;
                    vector<int> allele_prev;
                    ref.resize(2);
                    allele_curr.resize(2);
                    allele_prev.resize(2);
                    for (int i = 0; i < 2; i++) {
                        ref[1-i] = r%n_references;
                        allele_curr[1-i] = haplotype_to_allele_curr[ref[1-i]];
                        allele_prev[1-i] = haplotype_to_allele_prev[ref[1-i]];
                        r = r / n_references;
                    }
                    if ((allele_curr[0] == -1) || (allele_curr[1] == -1)) {continue;}     // Skipping nodes where the allele is not defined.
                    unsigned int index_prev = current_indexer->get_index(b_index, r_index);     // Getting the index of the node from the column column_index. This is needed for the beta value update.
                    long double beta_helper_0;
                    if ((allele_prev[0] == -1) || (allele_prev[1] == -1)) {
                        beta_helper_0 = 0.0L;   // If the corresponding node in previous column has unknown allle, make the helper variable from that as 0.
                    }
                    else {
                        beta_helper_0 = previous_projection_column->at(index_prev) * emission_probability_computer.at(allele_prev[0],allele_prev[1]);     // Calculating the beta(v,R1,R2).em(A1,A2) value from the notes. Since this value is used multiple times, just calculating it once here as a helper variable. Here we have to access allele from previous column.
                    }
                    // Updating the beta value
                    current_projection_column->at(index+r_index) += pow(current_transition_table->qr,2)*beta_helper_0 + current_transition_table->qr*current_transition_table->pr*(beta_helper_2[ref[0]]+beta_helper_3[ref[1]]-(2*beta_helper_0)) + pow(current_transition_table->pr,2)*(beta_helper_1-beta_helper_2[ref[0]]-beta_helper_3[ref[1]]+beta_helper_0);
                }
            }
        }
    }
    else {
        // Code block to calculate the sum of previous_projection_column and normalise backward_pass_column_table[0]
        for (int index = 0; index < previous_projection_column->size(); index++) scaling_sum += previous_projection_column->at(index);
    }
    // go through (old) projection column to scale the values -> when we lookup betas later, they will sum up to 1
    if(previous_projection_column != 0){
        std::transform((*previous_projection_column).begin(), (*previous_projection_column).end(), (*previous_projection_column).begin(), std::bind2nd(std::divides<long double>(), scaling_sum));
    }
    if(current_projection_column != 0){
        std::transform((*current_projection_column).begin(), (*current_projection_column).end(), (*current_projection_column).begin(), std::bind2nd(std::divides<long double>(), scaling_sum));
        backward_pass_column_table[column_index-1] = current_projection_column;
    }
    scaling_parameters[column_index] = scaling_sum;
}

// given the current matrix column, compute the forward probability table
void GenotypeHMM::compute_forward_column(size_t column_index, unique_ptr<vector<const Entry*>> current_input_column)
{
    // IMPORTANT: Here we are calculating the forward_column at column_index (and not column_index - 1 like with backward pass)!
    assert(column_index < input_column_iterator.get_column_count());

    Column* current_indexer = hmm_columns[column_index];
    assert(current_indexer != nullptr);

    // if the current input column was not provided, then create it
    if(current_input_column.get() == nullptr) {
        input_column_iterator.jump_to_column(column_index);
        current_input_column = input_column_iterator.get_next();
    }

    // obtain previous projection column (which is assumed to have already been computed)
    vector<long double>* previous_projection_column = nullptr;
    if (column_index > 0) {
        // forward_pass_column_table is a vector (containing vectors) of size 1 (since we only need to store one column at a time). Jana plans on changing it to a vector containing long double values.
        previous_projection_column = forward_pass_column_table[0];
        assert(previous_projection_column != nullptr);
    }

    // obtain the backward projection table, from where to get the backward probabilities
    // ASSUMED THAT THIS WORKS
    size_t k = (size_t)sqrt(input_column_iterator.get_column_count());
    vector<long double>* backward_probabilities = nullptr;
    backward_probabilities = backward_pass_column_table[column_index];
    // if column is not stored, recompute it
    if(backward_probabilities == nullptr){
        // compute index of next column that has been stored
        size_t next = std::min((unsigned int) ( ((column_index + k) / k) * k ), input_column_iterator.get_column_count()-1);
        for(size_t i = next; i > column_index; --i){
            transition_probability_table[i-1] = new TransitionProbabilityComputer(recombcost[i-1], allele_references->at(i));
            compute_backward_column(i);
        }
        if (backward_pass_column_table[column_index] ==  nullptr) {
            compute_backward_column(next);
            delete backward_pass_column_table[next-1];
            backward_pass_column_table[next-1] = nullptr;
        }
        assert (backward_pass_column_table[column_index] != nullptr);
        // last column just computed still needs to be scaled
        std::transform((*backward_pass_column_table[column_index]).begin(), (*backward_pass_column_table[column_index]).end(), (*backward_pass_column_table[column_index]).begin(), std::bind2nd(std::divides<long double>(), scaling_parameters[column_index]));
    }
    backward_probabilities = backward_pass_column_table[column_index];
    assert(backward_probabilities != nullptr);
    // initialize the new projection column (2D: has entry for every bipartition and transmission value)
    vector<long double>* current_projection_column = nullptr;
    if(column_index < input_column_iterator.get_column_count()){
        current_projection_column = new vector<long double>(current_indexer->get_column_size(), 0.0L);
    }
    int n_alleles = variant_n_allele_positions->at(column_index);
    Vector2D<long double> emission_probability_computer = Vector2D<long double>(n_alleles, n_alleles);
    TransitionProbabilityComputer* current_transition_table;
    if (column_index > 0) {
        current_transition_table = transition_probability_table[column_index - 1];
    }
    // sum of alpha*beta, used to normalize the likelihoods
    long double normalization = 0.0L;
    long double scaling_sum = 0.0L;
    long double sum = 0.0L;
    int b_index;
    int allele_1;
    int allele_2;
    unsigned int r_index;
    unsigned int index;
    vector<unsigned int> compatible_bipartitions;
    vector<int> haplotype_to_allele_curr = allele_references->at(column_index);     // This contains the haplotype-to-allele mapping for the position column_index
    vector<int> haplotype_to_allele_prev;
    if (column_index > 0) {
        haplotype_to_allele_prev = allele_references->at(column_index-1);     // This contains the haplotype-to-allele mapping for the position column_index-1
    }
    // iterate over all bipartitions
    unique_ptr<ColumnIndexingIterator> iterator = current_indexer->get_iterator(read_set);
    while (iterator->has_next()) {
        int bit_changed = -1;
        iterator->advance(&bit_changed);
        // Update the emission probability based on the bipartition defined by the iterator
        update_emission_probability(&emission_probability_computer, bit_changed, *iterator, *current_input_column);
        /* if (column_index == 33) cout << emission_probability_computer << endl;
        cout << endl;
        if (column_index == 33) {
            for (auto e: *current_input_column) cout << *e << endl;
        } */
        // Determine the indices that are compatible with with the bipartition at position column_index
        b_index = iterator->get_b_index();
        // Calculating the current_projection_column
        if (column_index == 0) {
            for (int r_index = 0; r_index < pow(n_references,2); r_index++) {
                // Extract the ref haplotype paths and alleles from r_index
                int r = r_index;    // Making a copy of r_index to extract R1, R2.
                vector<unsigned int> ref;
                vector<int> allele;
                ref.resize(2);
                allele.resize(2);
                for (int i = 0; i < 2; i++) {
                    ref[1-i] = r%n_references;
                    allele[1-i] = haplotype_to_allele_curr[ref[1-i]];
                    r = r / n_references;
                }
                index = current_indexer->get_index(b_index, r_index);
                if ((allele[0] == -1) || (allele[1] == -1)) {
                    current_projection_column->at(index) = 0;     // Skipping nodes where the allele is not defined.
                    continue;
                }
                current_projection_column->at(index) = emission_probability_computer.at(allele[0], allele[1]);  // For the first column, the alpha value is just emission probability based.
            }
        }
        else {
            compatible_bipartitions = hmm_columns[column_index-1]->get_backward_compatible_bipartitions(b_index, read_set);
            // Looping through all the compatible bipartitions in column column_index - 1.
            // Here b is the bipartition index of the previous column!
            // The bipartition index of the current column is stored in b_index.
            for (unsigned int b : compatible_bipartitions) {
                long double alpha_helper_1 = 0.0L;       // This helper value is the alpha(*,*) value.
                vector<long double> alpha_helper_2(n_references, 0.0L);      // This helper value is the alpha(R1,*) value.
                vector<long double> alpha_helper_3(n_references, 0.0L);      // This helper value is the alpha(*,R2) value.
                // Looping through all the nodes in bipartition b of column_index - 1 to calculate helper variables.
                for (int r_index = 0; r_index < pow(n_references, 2); r_index++) {
                    // Extract the ref haplotype paths and alleles from r_index
                    int r = r_index;    // Making a copy of r_index to extract R1, R2.
                    vector<unsigned int> ref;
                    vector<int> allele;
                    ref.resize(2);
                    allele.resize(2);
                    for (int i = 0; i < 2; i++) {
                        ref[1-i] = r%n_references;
                        allele[1-i] = haplotype_to_allele_prev[ref[1-i]];
                        r = r / n_references;
                    }
                    if ((allele[0] == -1) || (allele[1] == -1)) {continue;}     // Skipping nodes where the allele is not defined.
                    index = (n_references*n_references*b)+r_index;
                    long double a = previous_projection_column->at(index);
                    // Updating the helper variables
                    alpha_helper_1 += a;
                    alpha_helper_2[ref[0]] += a;
                    alpha_helper_3[ref[1]] += a;
                }
                // Looping through all the nodes in bipartition b_index in column_index. Will add all the alpha values in the bipartition coming from bipartition b of column_index -1.
                unsigned int index_prev = b * pow(n_references,2);      // Base index for the previous bipartition. So this index corresponds to the (0,0) ref haplotype of the bipartition. So to access the rest of the indices, we need to sum index_prev and r_index.
                for (int r_index = 0; r_index < pow(n_references, 2); r_index++) {
                    // Extract the ref haplotype paths and alleles from r_index
                    int r = r_index;    // Making a copy of r_index to extract R1, R2.
                    vector<unsigned int> ref;
                    vector<int> allele_prev;
                    vector<int> allele_curr;
                    ref.resize(2);
                    allele_prev.resize(2);
                    allele_curr.resize(2);
                    for (int i = 0; i < 2; i++) {
                        ref[1-i] = r%n_references;
                        allele_prev[1-i] = haplotype_to_allele_prev[ref[1-i]];
                        allele_curr[1-i] = haplotype_to_allele_curr[ref[1-i]];
                        r = r / n_references;
                    }
                    if ((allele_curr[0] == -1) || (allele_curr[1] == -1)) {continue;}     // Skipping nodes where the allele in column column_index is not defined.
                    index = current_indexer->get_index(b_index, r_index); // Here I am using ->get_index() since I want to find the index in the current column. In the previous for loop, I did not use this since current_indexer is not defined for column_index - 1.
                    long double alpha_helper_0;
                    if ((allele_prev[0] == -1) || (allele_prev[1] == -1)) {
                        alpha_helper_0 = 0.0L;   // If the corresponding node in previous column has unknown allle, make the helper variable from that as 0.
                    }
                    else {
                        alpha_helper_0 = previous_projection_column->at(index_prev+r_index);        // Storing alpha value of the compatible biparition of column_index - 1 having the same r_index as the bipartition in column_index. (Strictly speaking, we dont need this variable but this makes the expression very similar to backward algorithm)
                    }
                    // Updating the alpha value
                    current_projection_column->at(index) += emission_probability_computer.at(allele_curr[0], allele_curr[1])*(pow(current_transition_table->qr,2)*alpha_helper_0 + current_transition_table->qr*current_transition_table->pr*(alpha_helper_2[ref[0]]+alpha_helper_3[ref[1]]-(2*alpha_helper_0)) + pow(current_transition_table->pr,2)*(alpha_helper_1-alpha_helper_2[ref[0]]-alpha_helper_3[ref[1]]+alpha_helper_0));
                }
            }
        }
    }
    long double forward_backward = 0.0L;
    
    assert (current_projection_column->size() == backward_probabilities->size());
    for (unsigned int i = 0; i < backward_probabilities->size(); i++) {
        vector<unsigned int> ref;
        ref.resize(2);
        unsigned int r_index = (unsigned int)(i%(int)pow(n_references,2));
        for (int j = 0; j < 2; j++) {
            ref[1-j] = r_index%n_references;
            r_index = r_index / n_references;
        }
        vector<unsigned int> alleles;
        alleles.push_back(allele_references->at(column_index).at(ref.at(0)));
        alleles.push_back(allele_references->at(column_index).at(ref.at(1)));
        // Get the genotype index
        unsigned int g_index = 0;
        for (int allele = 0; allele < 2; allele++) {
            g_index += binomial_coefficient(allele + alleles.at(allele), alleles.at(allele) - 1);
        }
        current_projection_column->at(i) = current_projection_column->at(i) / scaling_parameters[column_index];
        sum += current_projection_column->at(i);
        forward_backward = current_projection_column->at(i) * backward_probabilities->at(i);
        // if (column_index == 33) cout << "Forward: " << current_projection_column->at(i) << "\tBackward: " << backward_probabilities->at(i) << "\tG Index: " << g_index << endl;
        normalization += forward_backward;
        
        // HARDCODED FOR A PEDIGREE SIZE OF 1.
        genotype_likelihood_table.at(0, column_index).likelihoods[g_index] += forward_backward;
    }
    /* if (column_index == 33) {
        cout << endl;
        for (auto g: genotype_likelihood_table.at(0, column_index).likelihoods) {
            cout << g << "\t";
        }
        cout << endl;
        exit(0);
    } */
    std::transform((*current_projection_column).begin(), (*current_projection_column).end(), (*current_projection_column).begin(), std::bind2nd(std::divides<long double>(), sum));
    // store the computed projection column (in case there is one)
    if(current_projection_column != 0){
        delete forward_pass_column_table[0];
        forward_pass_column_table[0] = current_projection_column;
    }
    // we can remove the backward-probability column
    if(backward_pass_column_table[column_index] != nullptr){
        delete backward_pass_column_table[column_index];
        backward_pass_column_table[column_index] = nullptr;
    }
    // scale the likelihoods
    for(size_t individuals_index = 0; individuals_index < pedigree->size(); ++individuals_index){
        genotype_likelihood_table.at(individuals_index,column_index).divide_likelihoods_by(normalization);
    }
}

vector<long double> GenotypeHMM::get_genotype_likelihoods(unsigned int individual_id, unsigned int position)
{
    assert(pedigree->id_to_index(individual_id) < genotype_likelihood_table.get_size0());
    assert(position < input_column_iterator.get_column_count());

    return genotype_likelihood_table.at(pedigree->id_to_index(individual_id),position).likelihoods;

}

void GenotypeHMM::update_emission_probability(Vector2D<long double>* em_prob, const int bit_changed, const ColumnIndexingIterator& iterator, vector<const Entry *>& entries) {
    int n_alleles = em_prob->get_size0();
    vector<int> binaryVector = iterator.get_binary_vector();
    assert (binaryVector.size() == entries.size());     // Asserting that the list of entries from current_input_column and the binary vector from current_indexer have the same size (Both are defined on the list of active reads).
    if (bit_changed >= 0) {
        // This is is executed when bipartitions are iterated and one bit is flipped.
        int newBit = binaryVector[bit_changed];     // Finding the new bit at the flip position
        if (entries.at(bit_changed)->get_allele_type() == -1) {
            // Not executing anything else if the allele type detected is unknown (this does not happen with the GAFs)
            return;
        }
        if (newBit == 1) {
            for (int i = 0; i < n_alleles; i++) {
                for (int j = 0; j < n_alleles; j++) {
                    // If the new_bit is 1, then the new value is such that the flipped read (entries.at(bit_changed)) has allele j (multiplication) and not allele i (division)
                    long double new_value = em_prob->at(i, j) * (entries.at(bit_changed)->get_emission_score()[j]/entries.at(bit_changed)->get_emission_score()[i]);
                    em_prob->set(i, j, new_value);}
            }
        }
        else {
            for (int i = 0; i < n_alleles; i++) {
                for (int j = 0; j < n_alleles; j++) {
                    // If the new_bit is 1, then the new value is such that the flipped read (entries.at(bit_changed)) has allele i (multiplication) and not allele j (division)
                    long double new_value = em_prob->at(i, j) * (entries.at(bit_changed)->get_emission_score()[i]/entries.at(bit_changed)->get_emission_score()[j]);
                    em_prob->set(i, j, new_value);
                }
            }
        }
    }
    else {
        // This is executed for the initialization step.
        for (int i = 0; i < n_alleles; i++) {
            for (int j = 0; j < n_alleles; j++) {
                // Filling the cell which indicates the bipartition 1 comes from allele i and bipartition 2 comes from allele j
                long double value = 1.0L;
                for (int entry_index = 0; entry_index < entries.size(); entry_index++) {
                    if (entries.at(entry_index)->get_allele_type() == -1) {
                        continue;
                    }
                    if (binaryVector[entry_index] == 0){
                        // If read at entry_index is in biparition 0 then value gets multiplied with emission from allele i
                        value = value * (entries.at(entry_index)->get_emission_score()[i]);
                    }
                    else {
                        assert (binaryVector[entry_index] == 1);
                        // If read at entry_index is in biparition 1 then value gets multiplied with emission from allele j
                        value = value * (entries.at(entry_index)->get_emission_score()[j]);
                    }
                }
                em_prob->set(i, j, value);
            }
        }
    }
}