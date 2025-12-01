/*
Taken from whatshap (version 2.8)
Original filename: src/columnindexingscheme.h
*/

#ifndef COLUMN_INDEXING_SCHEME_H
#define COLUMN_INDEXING_SCHEME_H

#include <vector>
#include <memory>
#include "phasingcolumnindexingiterator.h"

class PhasingColumnIndexingIterator;

class PhasingColumnIndexingScheme {
private:
	std::vector<uint32_t> read_ids;
	const PhasingColumnIndexingScheme* previous_column;
	const PhasingColumnIndexingScheme* next_column;
	uint32_t backward_projection_width;
	uint32_t forward_projection_width;
	std::vector<uint32_t>* forward_projection_mask;

public:
	
	/** Constructor.
	 * @param previousReadIDs IDs of reads active
	 */
	PhasingColumnIndexingScheme(const PhasingColumnIndexingScheme* previous_column, const std::vector<uint32_t>& read_ids);

	~PhasingColumnIndexingScheme();

	/** Set pointer to indexing scheme of next column. MUST be called before the getIterator
	 *  method is called. */
	void set_next_column(const PhasingColumnIndexingScheme* next_column);

	std::unique_ptr<PhasingColumnIndexingIterator> get_iterator();

	uint32_t column_size();

	uint32_t forward_projection_size();
 
	uint32_t get_forward_projection_width();

	uint32_t get_backward_projection_width();

	// return a const pointer to the read ids
	const std::vector<uint32_t> * get_read_ids();

	// return const forward projection mask (for debugging)
	const std::vector<uint32_t> * get_forward_projection_mask();

	friend class PhasingColumnIndexingIterator;
};

#endif
