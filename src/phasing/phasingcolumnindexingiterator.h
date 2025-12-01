/*
Taken from whatshap (version 2.8)
Original filename: src/columnindexingiterator.h
*/

#ifndef PHASING_COLUMN_INDEXING_ITERATOR_H
#define PHASING_COLUMN_INDEXING_ITERATOR_H

#include "../graycodes.h"

class PhasingColumnIndexingScheme;

class PhasingColumnIndexingIterator {
	
	private:
		const PhasingColumnIndexingScheme* parent;
		GrayCodes* graycodes;
		uint32_t index;
		uint32_t forward_projection;

	public:
		PhasingColumnIndexingIterator(const PhasingColumnIndexingScheme* parent);
		virtual ~PhasingColumnIndexingIterator();

		bool has_next();

		/** Move to next index (i.e. DP table row).
		 *
		 *  @param bit_changed If not null, and only one bit in the
		 *  partitioning (as retrieved by get_partition) is changed by this
		 *  call to advance, then the index of this bit is written to the
		 *  referenced variable; if not, -1 is written.
		 */
		void advance(int* bit_changed = 0);

		/** Index of the projection of the current read set onto the intersection between current and next read set. */
		uint32_t get_forward_projection();

		/** Index of the projection of the current read set onto the intersection between previous and the current read set. */
		uint32_t get_backward_projection();

		/** Row index in the DP table (within the current column). */
		uint32_t get_index();

		/** Bit-wise representation of the partitioning corresponding to the current index. */
		uint32_t get_partition();

		/** get index's backward projection (given index i), so that we don't have to iterate up to it, just to get it */
		uint32_t index_backward_projection(uint32_t i);

		/** get index's forward projection */
		uint32_t index_forward_projection(uint32_t i);
};

#endif
