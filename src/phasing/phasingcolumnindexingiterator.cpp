/*
Taken from whatshap (version 2.8)
Original filename: src/columnindexingiterator.cpp
*/

// DONE

#include <cassert>
#include "phasingcolumnindexingscheme.h"
#include "phasingcolumnindexingiterator.h"

using namespace std;

PhasingColumnIndexingIterator::PhasingColumnIndexingIterator(const PhasingColumnIndexingScheme* parent) {
	assert(parent != 0);
	this->parent = parent;
	this->graycodes = new GrayCodes(parent->read_ids.size());
	this->index = -1;
	this->forward_projection = -1;
}


PhasingColumnIndexingIterator::~PhasingColumnIndexingIterator() {
	delete graycodes;
}


bool PhasingColumnIndexingIterator::has_next() {
	return graycodes->has_next();
}


void PhasingColumnIndexingIterator::advance(int* bit_changed) {
	assert(graycodes->has_next());

	int graycode_bit_changed = -1;
	index = graycodes->get_next(&graycode_bit_changed);
	// first iteration?
	if (graycode_bit_changed == -1) {
		assert(index == 0);
		if (parent->forward_projection_mask != 0) {
			forward_projection = 0;
		}
	} else {
		if (parent->forward_projection_mask != 0) {
			// index of bit in the forward_projection
			int bit_index = parent->forward_projection_mask->at(graycode_bit_changed);
			if (bit_index >= 0) {
				forward_projection = forward_projection ^ (((uint32_t)1) << bit_index);
			}
		}
	}
	if (bit_changed != 0) {
		*bit_changed = graycode_bit_changed;
	}
}


uint32_t PhasingColumnIndexingIterator::get_forward_projection() {
	assert(index >= 0);
	return forward_projection;
}


uint32_t PhasingColumnIndexingIterator::get_backward_projection() {
	assert(index >= 0);
	return index & ((((uint32_t)1)<<parent->backward_projection_width) - 1);
}


uint32_t PhasingColumnIndexingIterator::get_index() {
	assert(index >= 0);
	return index;
}


uint32_t PhasingColumnIndexingIterator::get_partition() {
	assert(index >= 0);
	return index;
}


uint32_t PhasingColumnIndexingIterator::index_backward_projection(uint32_t i) {
	assert(i >= 0); // assert the proper boundaries
	assert(i < (((uint32_t)1) << parent->read_ids.size()));

	return i & ((((uint32_t)1) << parent->backward_projection_width) -1);
}


uint32_t PhasingColumnIndexingIterator::index_forward_projection(uint32_t i) {
	assert(i >= 0);
	assert(i < (((uint32_t)1) << parent->read_ids.size()));

	uint32_t i_forward_projection = 0;
	uint32_t s = 1;
	for(int j=0; j< parent->read_ids.size(); ++j) {
		uint32_t m = parent->forward_projection_mask->at(j);
		if(m != -1) {
			uint32_t s = (((uint32_t)1) << m);
			i_forward_projection += (s&i);
		}
	}

	return i_forward_projection;
}
