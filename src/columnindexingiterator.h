// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

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
		uint32_t b_index;
		/**
		 * The positions of the read clusters in the binary vector which are free to vary in Gray Code ordering.
		 */
		std::vector<uint32_t> free_positions;
		/**
		 * Constrained position mapping
		 * Key is the position of the representative read cluster in the binary vector (given by the min of the constrained pair)
		 * Value is the position of the read cluster which is constrained with the representative in the binary vector (max of the constrained pair)
		 */
		std::unordered_map<uint32_t, uint32_t> constrained_position_map;
		
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

		/**
		 * Check if the bit changed corresponds to a cluster
		 */
		bool is_clustered_bit(uint32_t bit_changed);

		/**
		 * Check if the bit changed corresponds to a read cluster that has a contrained cluster.
		 */
		bool has_constrained_bit(uint32_t bit_changed);

		/**
		 * For the representative constrained position given by bit_changed,
		 * we find the non-representative cluster position.
		 */
		uint32_t get_constrained_bit(uint32_t bit_changed);

		/**
		 * Get read ids of a cluster at the position
		 */
		std::vector<uint32_t> get_reads_from_cluster_id(uint32_t cluster_id);

		// Returns the index for the current bipartition of untagged active reads IDs that the iterator has processed
		uint32_t get_b_index();

		/**
		 * For a bit changed during the Gray Code ordering, if clusters are involved, the multiple reads have their bits flipped.
		 * Each cluster will have multiple reads and all of its bipartition gets flipped.
		 * If there are constrained clusters, then the reads for the constrained cluster will be flipped as well.
		 * 
		 * The input but is the bit in the Gray Code ordering that was changed.
		 * 
		 * The output is an unordered map of
		 *  - Key: read_id that has been flipped.
		 *  - Value: the new bit
		 * 
		 * This function generalises the cases of unclustered reads and clusters which do not have constraints.
		 */
		std::unordered_map<uint32_t, bool> get_changed_bits(uint32_t bit_changed);
	};

#endif
