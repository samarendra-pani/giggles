/*
Taken from whatshap (version 2.8)
Original filename: src/columniterator.cpp
*/

#include <cassert>
#include <limits>
#include <unordered_set>
#include <stdexcept>

#include "phasingcolumniterator.h"

using namespace std;

PhasingColumnIterator::PhasingColumnIterator(const ReadSet& set, const std::vector<variant_information_t>* variant_info_table, bool is_first_phasing_round) : 
	set(set), 
	is_first_phasing_round(is_first_phasing_round),
	variant_info_table(variant_info_table) {
	
	n = 0;
	next_read_index = 0;
	// create a mapping of genomic positions to column indices
	std::unordered_map<uint32_t, size_t> position_map;
	for (size_t i=0; i<variant_info_table->size(); ++i) {
		position_map[variant_info_table->at(i).position] = i;
	}
	// precompute first_reads
	first_reads.assign(variant_info_table->size(),  numeric_limits<size_t>::max());
	int pos = 0;
	for (size_t i=0; i<set.size(); ++i) {
		const Read* read = set.get(i);
		if (read->firstPosition() < pos) {
			throw std::runtime_error("PhasingColumnIterator: reads in ReadSet are not sorted.");
		}
		if (!read->isSorted()) {
			throw std::runtime_error("PhasingColumnIterator: encountered read with unsorted variants.");
		}
		/**
		 * @note IS THIS CORRECT?
		 */
		if (!read->isSelected()) {
			// Skipping unselected reads.
			continue;
		}
		auto first_column_it = position_map.find(read->firstPosition());
		auto last_column_it = position_map.find(read->lastPosition());
		assert(first_column_it != position_map.end());
		assert(last_column_it != position_map.end());
		assert(first_column_it->second <= last_column_it->second);
		assert(last_column_it->second < variant_info_table->size());
		for (size_t j=first_column_it->second; j<=last_column_it->second; ++j) {
			if (first_reads[j] == numeric_limits<size_t>::max()) {
				first_reads[j] = i;
			}
		}
		pos = read->firstPosition();
	}
	// For positions not covered by any read, fill in the index of next read that will
	// become active in subsequent columns
	if (first_reads.size() >= 2) {
		size_t next_index = first_reads[first_reads.size()-1];
		for (int i=first_reads.size()-2; i>=0; --i) {
			if (first_reads[i] == numeric_limits<size_t>::max()) {
				first_reads[i] = next_index;
			} else {
				next_index = first_reads[i];
			}
		}
	}
}


PhasingColumnIterator::~PhasingColumnIterator() {
	for (size_t i=0; i<blank_entries.size(); ++i) {
		delete blank_entries[i];
	}
	blank_entries.clear();
}


uint32_t PhasingColumnIterator::get_column_count() {
	return variant_info_table->size();
}


uint32_t PhasingColumnIterator::get_read_count() {
	return set.size();
}


const uint32_t PhasingColumnIterator::get_position(uint32_t i) {
	assert(i < variant_info_table->size());
	return variant_info_table->at(i).position;
}


bool PhasingColumnIterator::has_next() {
	return n < variant_info_table->size();
}


unique_ptr<vector<const Entry*> > PhasingColumnIterator::get_next() {
	// genomic position of the column to be returned
	int next_pos = variant_info_table->at(n).position;
	// check which of the current reads remain active
	list<active_read_t>::iterator list_it = active_reads.begin();
	while (list_it != active_reads.end()) {
		const Read* read = set.get(list_it->read_index);
		if (read->lastPosition() < next_pos) {
			list_it = active_reads.erase(list_it);
			continue;
		}
		while (read->getPosition(list_it->active_entry) < next_pos) {
			list_it->active_entry += 1;
			assert(list_it->active_entry < read->getVariantCount());
		}
		++list_it;
	}

	// check which new reads might become active
	while (next_read_index < set.size()) {
		const Read* read = set.get(next_read_index);
		if (read->isSelected() == false) {
			// skip unselected reads
			next_read_index += 1;
			continue;
		}
		int read_start = read->firstPosition();
		if (read_start == next_pos) {
			active_reads.push_back(active_read_t(next_read_index));
			next_read_index += 1;
		} else {
			assert(read_start > next_pos);
			break;
		}
	}

	// gather entries from active reads
	unique_ptr<vector<const Entry*> > result(new vector<const Entry*>());
	for (list_it = active_reads.begin(); list_it != active_reads.end(); ++list_it) {
		Read* read = set.get(list_it->read_index);
		if (!variant_info_table->at(n).phasable) {
			// the position has multiple possible alleles.
			// cannot phase
			Entry* e = new Entry(read->getID(), std::vector<uint32_t>{});
			blank_entries.push_back(e);
			result->push_back(e);
			continue;
		}
		if (is_first_phasing_round && variant_info_table->at(n).is_sv) {
			// in the first phasing round, we do not consider structural variants
			Entry* e = new Entry(read->getID(), std::vector<uint32_t>{});
			blank_entries.push_back(e);
			result->push_back(e);
			continue;
		}
		// Does read cover the current position?
		if (read->getPosition(list_it->active_entry) == next_pos) {
			// If so, add the entry to the result is the entry is biallelic
			Entry* entry = read->getEntry(list_it->active_entry);
			if (!entry->has_allele_type()) { entry->set_allele_type(variant_info_table->at(n).active_alleles); }
			result->push_back(entry);
		} else {
			// if not, generate a blank entry
			Entry* e = new Entry(read->getID(), std::vector<uint32_t>{});
			blank_entries.push_back(e);
			result->push_back(e);
		}
	}

	n += 1;
	return result;
}


void PhasingColumnIterator::jump_to_column(size_t k) {
	if (k == n) return;
	assert(k < variant_info_table->size());
	active_reads.clear();
	n = k;
	next_read_index = first_reads[k];
	u_int32_t pos = variant_info_table->at(k).position;

	// determine set of active reads
	while (next_read_index < set.size()) {
		const Read* read = set.get(next_read_index);
		if (read->lastPosition() < pos) {
			next_read_index += 1;
			continue;
		}
		if (read->firstPosition() <= pos) {
			size_t active_entry = 0;
			while (read->getPosition(active_entry) < pos) {
				active_entry += 1;
				assert(active_entry < read->getVariantCount());
			}
			active_reads.push_back(active_read_t(next_read_index, active_entry));
			next_read_index += 1;
		} else {
			break;
		}
	}
}
