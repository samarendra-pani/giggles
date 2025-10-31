/*
Taken from whatshap (version 2.8)
Original filename: src/columnindexingscheme.cpp
*/

// DONE

#include <cassert>
#include "phasingcolumnindexingiterator.h"
#include "phasingcolumnindexingscheme.h"

using namespace std;

PhasingColumnIndexingScheme::PhasingColumnIndexingScheme(const PhasingColumnIndexingScheme* previous_column, const std::vector<unsigned int>& read_ids) : read_ids(read_ids) {
	this->previous_column = previous_column;
	this->next_column = 0;
	// assert that read ids are ordered
	if (read_ids.size() > 0) {
		for (size_t i=0; i<read_ids.size()-1; ++i) {
			assert(read_ids[i] < read_ids[i+1]);
		}
	}
	this->forward_projection_mask = 0;
	this->backward_projection_width = 0;
	this->forward_projection_width = 0;
	if (previous_column != 0) {
		int i = 0;
		int j = 0;
		while ((i<previous_column->read_ids.size()) && (j<read_ids.size())) {
			if (previous_column->read_ids[i] == read_ids[j]) {
				backward_projection_width += 1;
				i += 1;
				j += 1;
			} else if (previous_column->read_ids[i] < read_ids[j]) {
				i += 1;
			} else {
				j += 1;
			}
		}
	}
}


PhasingColumnIndexingScheme::~PhasingColumnIndexingScheme() {
	if (forward_projection_mask != 0) delete forward_projection_mask;
}


unsigned int PhasingColumnIndexingScheme::column_size() {
	return ((unsigned int)1) << read_ids.size();
}


unsigned int PhasingColumnIndexingScheme::forward_projection_size() {
	return ((unsigned int)1) << forward_projection_mask->size();
}


unsigned int PhasingColumnIndexingScheme::get_forward_projection_width() {
	return forward_projection_width;
}


unsigned int PhasingColumnIndexingScheme::get_backward_projection_width() {
	return backward_projection_width;
}


void PhasingColumnIndexingScheme::set_next_column(const PhasingColumnIndexingScheme* next_column) {
	assert(next_column != 0);
	this->next_column = next_column;
	if (forward_projection_mask != 0) delete forward_projection_mask;
	forward_projection_width = 0;
	forward_projection_mask = new vector<unsigned int>(read_ids.size(),-1);
	int i = 0;
	int j = 0;
	int n = 0;
	while ((i<next_column->read_ids.size()) && (j<read_ids.size())) {
		if (next_column->read_ids[i] == read_ids[j]) {
			forward_projection_mask->at(j) = n;
			n += 1;
			i += 1;
			j += 1;
		} else if (next_column->read_ids[i] < read_ids[j]) {
			i += 1;
		} else {
			j += 1;
		}
	}

	forward_projection_width = n+1;
}


unique_ptr<PhasingColumnIndexingIterator> PhasingColumnIndexingScheme::get_iterator() {
	return unique_ptr<PhasingColumnIndexingIterator>(new PhasingColumnIndexingIterator(this));
}


const vector<unsigned int> * PhasingColumnIndexingScheme::get_read_ids() {
	return &(this->read_ids);
}


const vector<unsigned int> * PhasingColumnIndexingScheme::get_forward_projection_mask() {
	return forward_projection_mask;
}
