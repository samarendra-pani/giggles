/*
Taken from whatshap (version 2.8)
Original filename: src/columnindexingscheme.h
*/

// DONE

#ifndef COLUMN_INDEXING_SCHEME_H
#define COLUMN_INDEXING_SCHEME_H

#include <vector>
#include <memory>
#include "phasingcolumnindexingiterator.h"

class PhasingColumnIndexingIterator;

class PhasingColumnIndexingScheme {
private:
	std::vector<unsigned int> read_ids;
	const PhasingColumnIndexingScheme* previous_column;
	const PhasingColumnIndexingScheme* next_column;
	unsigned int backward_projection_width;
	unsigned int forward_projection_width;
	std::vector<unsigned int>* forward_projection_mask;

public:
	
	/** Constructor.
	 * @param previousReadIDs IDs of reads active
	 */
	PhasingColumnIndexingScheme(const PhasingColumnIndexingScheme* previous_column, const std::vector<unsigned int>& read_ids);

	~PhasingColumnIndexingScheme();

	/** Set pointer to indexing scheme of next column. MUST be called before the getIterator
	 *  method is called. */
	void set_next_column(const PhasingColumnIndexingScheme* next_column);

	std::unique_ptr<PhasingColumnIndexingIterator> get_iterator();

	unsigned int column_size();

	unsigned int forward_projection_size();
 
	unsigned int get_forward_projection_width();

	unsigned int get_backward_projection_width();

	// return a const pointer to the read ids
	const std::vector<unsigned int> * get_read_ids();

	// return const forward projection mask (for debugging)
	const std::vector<unsigned int> * get_forward_projection_mask();

	friend class PhasingColumnIndexingIterator;
};

#endif
