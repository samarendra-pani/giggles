// Code taken from WhatsHap (https://github.com/whatshap/whatshap)

#include <cassert>
#include <limits>
#include <unordered_set>
#include <stdexcept>

#include "columniterator.h"

using namespace std;

ColumnIterator::ColumnIterator(const ReadSet& set, const std::vector<variant_information_t>* variant_info_table) : set(set), variant_info_table(variant_info_table) {
	if(variant_info_table->size() == 0) return;
	this->m = variant_info_table->size()-1;
	this->n = 0;
	this->next_read_index = 0;

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
			throw std::runtime_error("ColumnIterator: reads in ReadSet are not sorted.");
		}
		if (!read->isSorted()) {
			throw std::runtime_error("ColumnIterator: encountered read with unsorted variants.");
		}
		auto first_column_it = position_map.find(read->firstPosition());
		auto last_column_it = position_map.find(read->lastPosition());
		assert(first_column_it != position_map.end());
		assert(last_column_it != position_map.end());
		assert(first_column_it->second <= last_column_it->second);
		assert(last_column_it->second <= variant_info_table->size());
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


ColumnIterator::~ColumnIterator() {
	for (size_t i=0; i<current_blank_entries.size(); ++i) {
		delete current_blank_entries[i];
	}
	current_blank_entries.clear();
}


uint32_t ColumnIterator::get_column_count() {
	return variant_info_table->size();
}


uint32_t ColumnIterator::get_read_count() {
	return set.size();
}


bool ColumnIterator::has_next() {
	return n < variant_info_table->size();
}

bool ColumnIterator::has_prev() {
	return m != (uint32_t)-1;
}


unique_ptr<vector<const Entry*>> ColumnIterator::get_next() {
	
	// clearing blank entries from previous column
	for (Entry* e: current_blank_entries) {
		delete e;
	}
	current_blank_entries.clear();

	// genomic position of the column to be returned
	uint32_t next_pos = variant_info_table->at(n).position;

	// getting active reads from column index n
	get_active_reads(n);
	
	list<active_read_t>::iterator list_it;
	// gather entries from active reads
	unique_ptr<vector<const Entry*> > result(new vector<const Entry*>());
	for (list_it = active_reads.begin(); list_it != active_reads.end(); ++list_it) {
		Read* read = set.get(list_it->read_index);
		if (!read->isClustered()) { continue; }
		// Does read cover the current position?
		if (read->getPosition(list_it->active_entry) == next_pos) {
			result->push_back(read->getEntry(list_it->active_entry));
		} 
		else {
			// if not, generate a blank entry
			Entry* e = new Entry(read->getID());
			current_blank_entries.push_back(e);
			result->push_back(e);
		}
	}
	n += 1;
	return result;
}

unique_ptr<vector<const Entry*>> ColumnIterator::get_prev() {
	
	// clearing blank entries from previous column
	for (Entry* e: current_blank_entries) {
		delete e;
	}
	current_blank_entries.clear();

	// genomic position of the column to be returned
	uint32_t next_pos = variant_info_table->at(m).position;

	// getting active reads from column index m
	get_active_reads(m);
	
	list<active_read_t>::iterator list_it;
	// gather entries from active reads
	unique_ptr<vector<const Entry*> > result(new vector<const Entry*>());
	for (list_it = active_reads.begin(); list_it != active_reads.end(); ++list_it) {
		Read* read = set.get(list_it->read_index);
		if (!read->isClustered()) { continue; }
		// Does read cover the current position?
		if (read->getPosition(list_it->active_entry) == next_pos) {
			result->push_back(read->getEntry(list_it->active_entry));
		} 
		else {
			// if not, generate a blank entry
			Entry* e = new Entry(read->getID());
			current_blank_entries.push_back(e);
			result->push_back(e);
		}
	}
	m -= 1;
	return result;
}


void ColumnIterator::jump_to_column(uint32_t k) {
	
	assert(k < variant_info_table->size());
	active_reads.clear();
	n = k;
	m = k;
}

void ColumnIterator::get_active_reads(uint32_t k) {
	active_reads.clear();
	next_read_index = first_reads[k];
	uint32_t pos = variant_info_table->at(k).position;

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
