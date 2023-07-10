#ifndef COLUMN_INDEXING_ITERATOR_H
#define COLUMN_INDEXING_ITERATOR_H

#include "graycodes.h"
#include "column.h"
#include "readset.h"

class Column;

class ColumnIndexingIterator {
private:
	Column* parent;
	GrayCodes* graycodes;
	unsigned int r_index;
	unsigned int b_index;
	// This contains the bipartition information for all the reads (0/1/-1)
	std::vector<int> binaryVector;
	bool hasBipartition;
	// This contains the position of reads whose haplotypes are unknown
	std::vector<int> freePositions;
	
public:
	ColumnIndexingIterator(Column* parent, ReadSet* set);
	virtual ~ColumnIndexingIterator();

	bool has_next();

	/** Move to next index (i.e. DP table row).
	  *
	  *  @param bit_changed If not null, and only one bit in the
	  *  partitioning (as retrieved by get_partition) is changed by this
	  *  call to advance, then the index of this bit is written to the
	  *  referenced variable; if not, -1 is written.
	  */
	void advance(int* bit_changed = 0);

	unsigned int get_b_index();

	std::vector<int> get_binary_vector() const;
};

#endif
