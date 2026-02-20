#ifndef BIPARTITIONITERATOR_H
#define BIPARTITIONITERATOR_H

#include "graycodes.h"
#include "column.h"
#include "readset.h"

/**
 * @brief Iterates over the various biparitions that can be formed by the read clusters at some variant position
 *
 * Based on the read clusters defined by the Column class, BipartitionIterator
 * identifies any constraints between clusters (as is possible with phased clusters).
 * It calls the GrayCodes class to provide the efficient ordering of the biparitions
 * of these clusters.
 * 
 * @see Column
 * @see GrayCodes
 */
class BipartitionIterator {
	
	private:
		/**
		 * The Column object which contains the information about read clusters
		 * at this position.
		 */
		Column* parent;
		
		/**
		 * GrayCodes object to iterate over the bipartitions of selected read clusters
		 */
		GrayCodes* graycodes;
		
		/**
		 * This stores the bit representation of which read cluster belongs to which bipartition.
		 * The read clusters comes from Column* parent.
		 * 
		 * So if the clusters are {1, 2, 3, 4}, then read_cluster_bit_representation = 5 will represent
		 * 
		 * bits -> ......| 0 | 1 | 0 | 1 |
		 * clusters ->   | 4 | 3 | 2 | 1 | <- Cluster 1 assigned to least significant bit.
		 * 
		 * So 4 and 2 are in biparition 0 and 3 and 1 are in biparition 1.
		 */
		uint32_t read_cluster_bit_representation;

		/**
		 * This stores the index that the GrayCode order gives. So this stores the info
		 * about the free positions.
		 * 
		 * The difference from read_cluster_bit_representation lies in the read clusters
		 * that are constraint with another. Constraints are represented as a MAX ID -> MIN ID map.
		 * 
		 * The MIN is representated in the free_positions vector and the MAX ID is changed based
		 * on what MIN gets assigned in Gray Code order. So MAX is not considered in GrayCode.
		 * So we want to remove its contribution.
		 * 
		 * Example:
		 * clusters -> 14 | 9 | 8 | 6 | 5 | 2 | 1
		 * state -> ..  F | C | F | C | F | F | F --- Constraint/Free
		 * bits -> ...  0 | 0 | 1 | 1 | 0 | 0 | 0 -> read_cluster_bit_representation
		 * 
		 * We remove the constriant positions to get
		 * clusters -> 14 | 8 | 5 | 2 | 1
		 * bits -> ...  0 | 1 | 0 | 0 | 0 -> bipartition index
		 * 
		 * Reason to keep this separate from read_cluster_bit_representation:
		 *  - This index goes from 0 to 2^(number of free positions - 1) - 1
		 *  - read_cluster_bit_representation is more chaotic since constrained pairs cause multiple bit flips.
		 * 
		 * So this index is ideal for memory allocations and access
		 */
		uint32_t bipartition_index;
		
		
	public:
		BipartitionIterator(Column* parent, ReadSet* set);
		virtual ~BipartitionIterator();

		bool has_next() const;

		/** 
		 * Move to next biparition of read clusters.
		 *
		 * @param cluster_bit_changed this parameter will store which cluster index was changed.
		 * At this position, there will be a vector of read clusters and the cluster index will point
		 * to the read cluster in this vector that was flipped.
		 * 
		 * In the first call, since we stay at all (free) clusters belonging to bipartition 0 (i.e. 
		 * there is no flipping), this variable stores -1.
		 */
		void advance(int* cluster_bit_changed);

		/**
		 * Check if the bit changed corresponds to a cluster
		 */
		bool is_clustered_bit(uint32_t cluster_bit_changed) const;

		/**
		 * Check if the bit changed corresponds to a read cluster that has a contrained cluster.
		 */
		bool has_constrained_bit(uint32_t cluster_bit_changed) const;

		/**
		 * For the representative constrained position given by bit_changed,
		 * we find the non-representative cluster position.
		 */
		uint32_t get_constrained_bit(uint32_t cluster_bit_changed) const;

		/**
		 * Get indices of read ids that are present in the cluster.
		 * The index is with respect to the parent column's read_ids.
		 */
		const std::vector<uint32_t>* get_read_index_from_cluster_id(uint32_t cluster_id) const;

		/**
		 * Returns pointer to the parent column.
		 */
		Column* get_parent_column() const;

		/**
		 * Return the bipartition we are currently working on.
		 * This biparition is formed using the Gray Code ordering of free positions.
		 */
		uint32_t get_bipartition_index() const;

		/**
		 * Returns the bit representation of the read clusters' bipartition.
		 */
		uint32_t get_read_cluster_bit_representation() const;

		/**
		 * For a bit changed during the Gray Code ordering, if clusters are involved, the multiple reads have their bits flipped.
		 * Each cluster will have multiple reads and all of its bipartition gets flipped.
		 * If there are constrained clusters, then the reads for the constrained cluster will be flipped as well.
		 * 
		 * The input but is the bit in the Gray Code ordering that was changed.
		 * 
		 * The output is an unordered map of
		 *  - Key: index of the read that has been flipped. The index is with respect to the parent column's read_ids.
		 *  - Value: the new bit
		 * 
		 * This function generalises the cases of unclustered reads and clusters which do not have constraints.
		 */
		void get_changed_reads(uint32_t cluster_bit_changed, std::unordered_map<uint32_t, bool>& changed_reads) const;
	};

#endif // BIPARTITIONITERATOR_H
