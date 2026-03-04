/*
Taken from whatshap (version 2.8)
Original filename: src/pedigreedptable.cpp
*/

#include <stdexcept>
#include <cassert>
#include <limits>
#include <fstream>
#include <array>
#include <algorithm>
#include <cmath>
#include <vector>

#include "phasingcolumncostcomputer.h"
#include "phasingdptable.h"
#include "readbipartitioning/phasesetcomputer.h"
#include "readbipartitioning/haplotagcomputer.h"
#include "readbipartitioning/set_cluster_ids.h"

using namespace std;

PhasingDPTable::PhasingDPTable(ReadSet* read_set, const vector<variant_information_t>* variant_info_table, bool first_phasing_round) :
	read_set(read_set),
	optimal_score(0u),
	optimal_score_index(0u),
	input_column_iterator(*read_set, variant_info_table, first_phasing_round),
	variant_info_table(variant_info_table)
{	
	compute_table();
	// creating the haplotypes as super reads
	ReadSet* superreads = new ReadSet();
	get_super_reads(superreads);
	// getting the optimal bipartition of reads used in the DP table
	const std::vector<bool> *optimal_partitioning = get_optimal_partitioning();
	// getting accessible positions
	std::vector<uint32_t>* accessible_positions = new std::vector<uint32_t>();
	for (uint32_t i = 0; i < variant_info_table->size(); ++i) {
		if (variant_info_table->at(i).count_active_alleles() <= 2) {
			assert(variant_info_table->at(i).count_active_alleles() == 2);
			accessible_positions->push_back(variant_info_table->at(i).position);
		}
	}
	compute_phasesets(accessible_positions, read_set, superreads);
	haplotag_selected_reads(read_set, optimal_partitioning);
	haplotag_unselected_reads(read_set, superreads); // phasesets have to be called before this function.
	set_read_cluster_ids(read_set);
	delete superreads;
	delete accessible_positions;
}


PhasingDPTable::~PhasingDPTable() {
	init(projection_column_table, 0);
	init(index_backtrace_table, 0);
	init(indexers, 0);
}


unique_ptr<vector<uint32_t> > PhasingDPTable::extract_read_ids(const vector<const Entry *>& entries) {
	unique_ptr<vector<uint32_t> > read_ids(new vector<uint32_t>());
	for (size_t i=0; i<entries.size(); ++i) {
		read_ids->push_back(entries[i]->get_read_id());
	}
	return read_ids;
}


size_t PhasingDPTable::popcount(size_t x) {
	uint32_t count = 0;
	for (;x; x >>= 1) {
		count += x & 1;
	}
	return count;
}


void PhasingDPTable::clear_table() {
	size_t column_count = input_column_iterator.get_column_count();

	init(projection_column_table, column_count);
	init(index_backtrace_table, column_count);
	init(indexers, column_count);

	index_path.clear();

	optimal_score = numeric_limits<uint32_t>::max();
	optimal_score_index = 0;
}


void PhasingDPTable::compute_table() {
	clear_table();

	// empty read-set, nothing to phase, so MEC score is 0
	if (input_column_iterator.get_column_count() == 0) {
		optimal_score = 0;
		optimal_score_index = 0;
		return;
	}

	input_column_iterator.jump_to_column(0);
	unique_ptr<vector<const Entry *> > current_input_column;
	unique_ptr<vector<const Entry *> > next_input_column;
	// get the next column ahead of time
	next_input_column = input_column_iterator.get_next();
	unique_ptr<vector<uint32_t> > next_read_ids = extract_read_ids(*next_input_column);
	PhasingColumnIndexingScheme* next_indexer = new PhasingColumnIndexingScheme(0, *next_read_ids);
	indexers[0] = next_indexer;

	// forward pass: create a sparse table, storing values at every sqrt(#columns)-th position,
	size_t k = (size_t)sqrt(input_column_iterator.get_column_count());
	for (size_t column_index=0; column_index<input_column_iterator.get_column_count(); ++column_index) {
		// make former next column the current one
		current_input_column = std::move(next_input_column);
		unique_ptr<vector<uint32_t> > current_read_ids = std::move(next_read_ids);

		PhasingColumnIndexingScheme* current_indexer = next_indexer;
		// peek ahead and get the next column
		if (input_column_iterator.has_next()) {
			next_input_column = input_column_iterator.get_next();
			next_read_ids = extract_read_ids(*next_input_column);
			next_indexer = new PhasingColumnIndexingScheme(current_indexer,*next_read_ids);
			current_indexer->set_next_column(next_indexer);
			indexers[column_index + 1] = next_indexer;
		} else {
			assert(next_input_column.get() == 0);
			assert(next_read_ids.get() == 0);
			next_indexer = 0;
		}

		compute_column(column_index, std::move(current_input_column));

		// determine whether to delete previous column (to save space)
		if ((k>1) && (column_index > 0) && (((column_index-1)%k) != 0)) {
			delete index_backtrace_table[column_index-1];
			delete projection_column_table[column_index-1];
			index_backtrace_table[column_index-1] = nullptr;
			projection_column_table[column_index-1] = nullptr;
		}
	}

	// perform a backtrace to get optimal path
	index_path.assign(indexers.size(), 0);
	uint32_t v;
	v = optimal_score_index;
	index_path[indexers.size()-1] = v;
	for(size_t i = indexers.size()-1; i > 0; --i) { // backtrack through table
		// ensure that index_backtrace_table[i-1] and transmission_backtrace_table[i-1] exist
		if (projection_column_table[i-1] == nullptr) {
			// compute index of last previous column that has been stored
			size_t j = (i-1) / k * k;
			assert(projection_column_table[j] != nullptr);
			for (j=j+1; j<i; ++j) {
				compute_column(j);
			}
		}
		// compute index value for the current column
		unique_ptr<PhasingColumnIndexingIterator> iterator = indexers[i]->get_iterator();
		uint32_t backtrace_index = iterator->index_backward_projection(v);
		v = index_backtrace_table[i-1]->at(backtrace_index);
		index_path[i-1] = v;
		// free parts of the DP table no longer needed
		if (i%k == 0) {
			for (size_t j=i; (j<i+k) && (j<input_column_iterator.get_column_count()-1); ++j) {
				assert(projection_column_table[j] != nullptr);
				delete index_backtrace_table[j];
				delete projection_column_table[j];
				index_backtrace_table[j] = nullptr;
				projection_column_table[j] = nullptr;
			}
		}
	}
}


void PhasingDPTable::compute_column(size_t column_index, unique_ptr<vector<const Entry*>> current_input_column) {
	assert(column_index < input_column_iterator.get_column_count());

	// check whether requested column is already there
	if (projection_column_table[column_index] != nullptr) {
		assert(index_backtrace_table[column_index] != nullptr);
		return;
	}

	PhasingColumnIndexingScheme* current_indexer = indexers[column_index];
	assert(current_indexer != nullptr);

	// if current input column was not provided, then create it
	if (current_input_column.get() == nullptr) {
		input_column_iterator.jump_to_column(column_index);
		current_input_column = input_column_iterator.get_next();
	}

	// reserve memory for the current DP column
	vector<uint32_t> dp_column(current_indexer->column_size(), 0);

	// obtain previous projection column (which is assumed to have been already computed)
	vector<uint32_t>* previous_projection_column = nullptr;
	if (column_index > 0) {
		previous_projection_column = projection_column_table[column_index - 1];
	}

	// initialize forward projection column and associated backtrace columns,
	// if existing (i.e. if not last column)
	vector<uint32_t>* current_projection_column = nullptr;
	vector<uint32_t>* index_backtrace_column = nullptr;
	if (column_index + 1 < input_column_iterator.get_column_count()) {
		current_projection_column = new vector<uint32_t>(
			current_indexer->forward_projection_size(),
			numeric_limits<uint32_t>::max()
		);
		index_backtrace_column = new vector<uint32_t>(
			current_indexer->forward_projection_size(),
			numeric_limits<uint32_t>::max()
		);
	}

	// create column cost computers
	PhasingColumnCostComputer cost_computer(*current_input_column, column_index, variant_info_table);
	
	// iterate over all bipartitions
	unique_ptr<PhasingColumnIndexingIterator> iterator = current_indexer->get_iterator();
	while (iterator->has_next()) {
		int bit_changed = -1;
		iterator->advance(&bit_changed);
		if (bit_changed >= 0) {
			cost_computer.update_partitioning(bit_changed);
		} else {
			cost_computer.set_partitioning(iterator->get_partition());
		}

		// Determine index in backward projection column from where to fetch the previous cost
		size_t backward_projection_index = 0;
		if (column_index > 0) {
			backward_projection_index = iterator->get_backward_projection();
		}
		// Determine index in the current DP column to be written
		size_t current_index = iterator->get_index();

		// Compute cost incurred by current cell of DP table
		uint32_t current_cost = cost_computer.get_cost();
		uint32_t min = numeric_limits<uint32_t>::max();
		// add up cost from current_cost column and previous columns
		uint32_t val;
		uint32_t previous_cost = 0;
		if (column_index > 0) {
			previous_cost = previous_projection_column->at(backward_projection_index);
		}
		if ((current_cost < numeric_limits<uint32_t>::max()) && (previous_cost < numeric_limits<uint32_t>::max())) {
			val = current_cost + previous_cost;
		} else {
			val = numeric_limits<uint32_t>::max();
		}

		// check for new minimum
		if (val < min) {
			min = val;
		}
		dp_column.at(current_index) = min;
		
		// if last DP column, then check for new optimal score, otherwise update forward projection and backtrace columns
		if (current_projection_column == 0) {
			// update running optimal score index
			if (dp_column.at(current_index) < optimal_score) {
				optimal_score = dp_column.at(current_index);
				optimal_score_index = iterator->get_index();
			}
		} else {
			uint32_t forward_index = iterator->get_forward_projection();
			uint32_t it_idx = iterator->get_index();
			if (dp_column.at(current_index) < current_projection_column->at(forward_index)) {
				current_projection_column->at(forward_index) = dp_column.at(current_index);
				index_backtrace_column->at(forward_index) = it_idx;
			}
		}
	}

	// if not last column, then store computed tables
	if (current_projection_column != 0) {
		index_backtrace_table[column_index] = index_backtrace_column;
		projection_column_table[column_index] = current_projection_column;
	}
}


uint32_t PhasingDPTable::get_optimal_score() {
	//if (backtrace_table.empty()) throw runtime_error("Empty backtrace table");
	return optimal_score;
}


void PhasingDPTable::get_super_reads(ReadSet* output_read_set) {
	assert(output_read_set != nullptr);
	assert(output_read_set->size() == 1);

	input_column_iterator.jump_to_column(0);
	
	std::pair<Read*,Read*> superreads;
	// removed sample id from the new Read declaration
	superreads = std::make_pair(
		new Read("superread_0", 0, (uint32_t)-1),
		new Read("superread_1", 0, (uint32_t)-1)
	);

	PhasingColumnCostComputer::phased_variant_t population_alleles;
	vector<uint32_t> active_alleles;
	uint32_t pos;
	uint32_t v;
	if (index_backtrace_table.empty()) {
		assert(!input_column_iterator.has_next());
	} else {
		// run through the file again with the input_column_iterator
		uint32_t i = 0; // column index
		while (input_column_iterator.has_next()) {
			v = index_path[i];
			unique_ptr<vector<const Entry *> > column = input_column_iterator.get_next();
			PhasingColumnCostComputer cost_computer(*column, i, variant_info_table);
			cost_computer.set_partitioning(v);

			population_alleles = cost_computer.get_alleles();
			// some sort of check if see if the alleles are blank?
			active_alleles = variant_info_table->at(i).get_active_positions();
			pos = input_column_iterator.get_position(i);
			if (active_alleles.size() > 2) {
				// This position was not phased. Adding BLANKs
				superreads.first->addVariant(pos, vector<long double>{});
				superreads.second->addVariant(pos, vector<long double>{});	
			} else {
				assert (active_alleles.size() == 2);
				// TODO: compute proper weights based on likelihoods.
				superreads.first->addVariant(pos, vector<uint32_t>(population_alleles.quality), population_alleles.allele0);
				superreads.second->addVariant(pos, vector<uint32_t>(population_alleles.quality), population_alleles.allele1);
			}
			++i; // next column
		}
	}
	assert(output_read_set != nullptr);
	output_read_set->add(superreads.first);
	output_read_set->add(superreads.second);
}


const vector<bool>* PhasingDPTable::get_optimal_partitioning() {
	vector<bool>* partitioning = new vector<bool>(read_set->size(),false);

	for(size_t i=0; i< index_path.size(); ++i) {
		uint32_t mask = 1; // mask to pass over the partitioning (i.e., index)
		for(size_t j=0; j< indexers[i]->get_read_ids()->size(); ++j) {
			uint32_t index = index_path[i];
			if((index & mask) == 0) { // id at this index is in p0 (i.e., in the part.)
				partitioning->at(indexers[i]->get_read_ids()->at(j)) = true;
			}
			mask = mask << 1;
		}
	}
	
	return partitioning;
}
