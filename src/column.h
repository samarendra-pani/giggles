// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#ifndef COLUMN_H
#define COLUMN_H

#include <vector>
#include <memory>
#include "columnindexingiterator.h"
#include "readset.h"

class ColumnIndexingIterator;

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
	 * Map cluster IDs to the reads that are present in this cluster at the variant position.
	 */
	std::unordered_map<uint32_t, std::vector<uint32_t>> cluster_id_to_read_ids_map;
	/**
	 * Read IDs that have been tagged during phasing are considered as "clusters" here.
	 * These read clusters are represented by a read ID defined during phasing.
	 * Untagged reads are represented by their own read IDs.
	 * Vector contains the clusters for the next variant position.
	 */
	std::vector<uint32_t> next_read_cluster_ids;
	/**
	 * Precomputed bipartitions using the read clusters.
	 * At runtime, just need the assignment of haplotype of the common clusters.
	 */
	std::vector<uint32_t> cached_bipartitions;
	/**
	 * Keepin track of the masks of the next read clusters based on the position on the cluster in the left column.
	 * maks = 1 << position_in_left_column
	 * By knowing these positions, we can quickly determine the compatible bipartitions between the two columns.
	 * For clusters that are unique to the next column, the mask will be 0.
	 */
	std::vector<uint32_t> next_read_cluster_masks;
	/**
	 * Keep track of the read cluster constraints (read clusters that cannot be in the same bipartition)
	 * This is stored as an unordered map from C_ID1 to C_ID2 (where the two clusters are constrained)
	 * Note that a mapping from C_ID2 to C_ID1 is not stored. So while doing Graycode ordering, we only shift the bits of C_ID1.
	 * We always store the mapping from max(C_ID1, C_ID2) to min(C_ID1, C_ID2).
	 * This is to avoid double counting of constraints.
	 */
	std::unordered_map<uint32_t, uint32_t> read_cluster_constraints;

	
public:

	Column(const uint32_t index, const std::vector<uint32_t>& read_ids, const std::vector<uint32_t>& next_read_ids, ReadSet* set);
	
	// return a pointer to the read cluster ids
	std::vector<uint32_t> * get_read_cluster_ids();

	// return a pointer to the read cluster constraints map
	std::unordered_map<uint32_t, uint32_t> * get_read_cluster_constraints_map();

	// return a pointer to the cluster id to read ids map
	std::unordered_map<uint32_t, std::vector<uint32_t>> * get_cluster_id_to_read_ids_map();

	// returns a pointer to the bipartition iterator which uses graycode.
	std::unique_ptr<ColumnIndexingIterator> get_iterator(ReadSet* set);

	// precompute bipartition compatibility information based on the read clusters
	void precompute_bipartition();

	// returns the compatible bipartitions of bipartition b_index (of pos v+1) in column v
	std::vector<uint32_t> get_backward_compatible_bipartitions(uint32_t b_index);

};

#endif
