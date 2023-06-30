#include <cassert>
#include "column.h"
#include "columnindexingiterator.h"
#include "readset.h"
#include "math.h"

using namespace std;

ColumnIndexingIterator::ColumnIndexingIterator(Column* parent, const ReadSet& set) {
	assert(parent != 0);
	this->parent = parent;
	
	int l = 0;
	for (int i = 0; i < parent->get_read_ids()->size(); i++) {
		binaryVector.push_back(set.get(i)->getHaplotag());
		if ( !set.get(parent->get_read_ids()->at(0))->hasHaplotag() ) {
			freePositions.push_back(i);
		}
	}
	// The Gray Code ordering is now done only on the positions whose bipartitions are unknown
	this->graycodes = new GrayCodes(freePositions.size());
	this->b_index = -1;
}


ColumnIndexingIterator::~ColumnIndexingIterator() {
	delete graycodes;
}


bool ColumnIndexingIterator::has_next() {
	return graycodes->has_next();
}


void ColumnIndexingIterator::advance(int* bit_changed) {
	assert(graycodes->has_next());
	int graycode_bit_changed = -1;
	unsigned int graycode_binaryindex = graycodes->get_next(&graycode_bit_changed);
	vector<int> graycode_binaryvector = graycodes->toBinary(graycode_binaryindex);
	
	// This indicates which bit in the gray code has been changed.
	// Since Gray Code is made on read subsets, this `graycode_bit_changed` has to be
	// translated to what bit has been changed in the entire set of active reads.
	if (bit_changed != 0) {
		*bit_changed = this->freePositions[graycode_bit_changed];
	}
	
	assert (graycode_binaryvector.size() == this->freePositions.size());

	// The main binary vector now has to be updated with the free positions.
	
	if (this->b_index == -1) {
		// Initially it is update element-wise and b-index accordingly calculated.
		this->b_index = 0;
		for (int i = 0; i < this->freePositions.size(); i++) {
			this->binaryVector[this->freePositions[i]] = graycode_binaryvector[i];
		}
		for (int i = 0; i < this->binaryVector.size(); i++) {
			this->b_index = this->b_index + this->binaryVector[i]*(pow(2,i));
		}
	}
	else {
		// For the subsequent updates to the bipartition index which can be easily done by single bit flips using masks.
		int new_bit = graycode_binaryvector[graycode_bit_changed];
		this-> binaryVector[*bit_changed] = new_bit;
		// Updating b_index using masks.
		int mask = 1 << *bit_changed;
		this->b_index ^= mask;
	}
}

unsigned int ColumnIndexingIterator::get_b_index() {
	return this->b_index;
}

vector<int> ColumnIndexingIterator::get_binary_vector() const {
	return this->binaryVector;
}
