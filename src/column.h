// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#ifndef COLUMN_H
#define COLUMN_H

#include <vector>
#include <memory>
#include "readset.h"

/**
 * Forward declaration of BipartitionIterator to avoid circular dependency between Column and BipartitionIterator.
 */
class BipartitionIterator;

/**
 * @brief Stores the information of active reads, clusters and auxilary data structures needed for connecting between Columns
 *
 * This class takes the information of which reads are active at some position and the next position.
 * It identifies all the read clusters and constraints between read clusters.
 * To pass to the BipartitionIterator, it finds all the clusters that can be iterated over and sorts
 *   them in a way to reduce number of operations for emission probability calculations.
 * 
 * It also stores some additional data structure for quickly creating connections between the next column and this one.
 * 
 * Cached Biparitions: These are created using the free read clusters that are unique to this column. When finding connections
 *   between the next column and this one, these read clusters can be in any bipartition. So we store all the gray code indices
 *   where the non-unique free read clusters are set to 0 and these unique ones are in all combinations.
 *   Later for quick retrieval of connections, we just need to set the bits of the non-unique free read clusters.
 * 
 * Next Read Cluster Masks: These are created using the free read clusters shared between the two columns (referred to as 
 *   non-unique free read clusters in "Cached Bipartitions"). Since we already know the bit positions of all these clusters
 *   in the gray code indices, we just to store these positions as masks and just multiply it whatever bipartition these clusters
 *   are in.
 * 
 * @see BipartitionIterator
 */
class Column {
	
	private:
		// read IDs of active reads at the variant position
		std::vector<uint32_t> read_ids;
		/**
		 * Read IDs that have been tagged during phasing are considered as "clusters" here.
		 * These read clusters are represented by a read ID defined during phasing.
		 * Untagged reads are represented by their own read IDs.
		 * Vector contains the clusters for variant position.
		 */
		std::vector<uint32_t> read_cluster_ids;
		/**
		 * Map cluster IDs to the read indices that are present in this cluster at the variant position.
		 */
		std::unordered_map<uint32_t, std::vector<uint32_t>> cluster_id_to_read_index_map;
		/**
		 * Precomputed bipartitions using the read clusters.
		 * At runtime, just need the assignment of haplotype of the common clusters.
		 */
		std::vector<uint32_t> cached_bipartitions;
		/**
		 * Keeping track of the masks of the next read clusters based on the position on the cluster in the left column.
		 * mask = 1 << (position of cluster in the gray code bit representation of left column)
		 * By knowing these positions, we can quickly determine the compatible bipartitions between the two columns.
		 * For clusters that are unique to the next column, the mask will be 0.
		 */
		std::vector<uint32_t> next_read_cluster_masks;
		/**
		 * The positions of the read clusters in the Gray Code index.
		 * The read clusters are sorted based on the number of reads associated with the cluster.
		 * The LSB has the least reads and the MSB has most reads.
		 */
		std::vector<uint32_t> sorted_free_read_cluster_positions;
		/**
		 * Constrained position mapping
		 * Key is the position of the representative read cluster (given by the min of the constrained pair)
		 * Value is the position of the non-representative read cluster (max of the constrained pair) in the read_cluster_bit_representation
		 * 
		 * This map is used to go from the representative to non-representative cluster to convert
		 * Gray Code index (which only has the representative cluster) to read_cluster_bit_representation (which has both clusters)
		 */
		std::unordered_map<uint32_t, uint32_t> constrained_position_map;

		// precompute bipartition compatibility information based on the read clusters
		void precompute_bipartition(std::vector<uint32_t>& next_read_cluster_ids);
		
	public:

		Column(const std::vector<uint32_t>& read_ids, const std::vector<uint32_t>& next_read_ids, ReadSet* set);
		
		// return a pointer to the read cluster ids
		const std::vector<uint32_t> * get_read_cluster_ids() const;

		// return a pointer to the read ids
		const std::vector<uint32_t> * get_read_ids() const;

		// return a pointer to the cluster id to read indices map
		const std::unordered_map<uint32_t, std::vector<uint32_t>> * get_cluster_id_to_read_index_map() const;

		// returns a pointer to the bipartition iterator which uses graycode.
		std::unique_ptr<BipartitionIterator> get_iterator(ReadSet* set);

		// returns the compatible bipartitions of bipartition read_cluster_bit_representation (of pos v+1) in column v
		void get_backward_compatible_bipartitions(uint32_t read_cluster_bit_representation, std::vector<uint32_t>& result) const;

		// returns the number of bipartitions this column has.
		uint32_t get_num_bipartition() const;

		// returns a pointer to the free read cluster positions
		const std::vector<uint32_t> * get_sorted_free_read_cluster_positions() const;

		// returns a pointer to the constrained position map
		const std::unordered_map<uint32_t, uint32_t> * get_constrained_position_map() const;

};

#endif
