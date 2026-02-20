// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#include <cassert>
#include <math.h>
#include <algorithm>
#include <numeric>

#include "column.h"

using namespace std;

Column::Column(const std::vector<uint32_t>& read_ids, const std::vector<uint32_t>& next_read_ids, ReadSet* set): 
	read_ids(read_ids) {
	
	std::vector<uint32_t> next_read_cluster_ids;
	/**
	 * Finding the read clusters from the phasing (and possibly other clustering later)
	 */
	uint32_t count = 0;
	Read* read_obj;
	uint32_t c_id;
	for (const auto& current_read_id : read_ids) {
        read_obj = set->get(current_read_id);

        if (!read_obj->isClustered()) {
			// Read is not clustered. So we use it as its own cluster.
            read_cluster_ids.push_back(current_read_id);
        } else {
            c_id = read_obj->getClusterID();
            
            // Only add unique cluster IDs to the list
            if (cluster_id_to_read_index_map.count(c_id) == 0) {
                read_cluster_ids.push_back(c_id);
            }
            cluster_id_to_read_index_map[c_id].push_back(count);
        }
		count ++;
    }
	/**
	 * Sorting the read clusters since they need not be pre-sorted.
	 */
	std::sort(read_cluster_ids.begin(), read_cluster_ids.end());
	/**
	 * Keep track of the read cluster constraints (read clusters that cannot be in the same bipartition)
	 * This is stored as an unordered map from C_ID1 to C_ID2 (where the two clusters are constrained)
	 * Note that a mapping from C_ID2 to C_ID1 is not stored. So while doing Graycode ordering, we only shift the bits of C_ID1.
	 * We always store the mapping from max(C_ID1, C_ID2) to min(C_ID1, C_ID2).
	 * This is to avoid double counting of constraints.
	 */
	std::unordered_map<uint32_t, uint32_t> read_cluster_constraints;
	uint32_t con_c_id;
	uint32_t first_read_idx;
	uint32_t dependent;
	uint32_t anchor;
	for (uint32_t c_id : read_cluster_ids) {
        // We can check just the first read of the cluster to find constraints
        // (Assuming all reads in a cluster share the constraint)
    	first_read_idx = cluster_id_to_read_index_map[c_id][0];
        read_obj = set->get(read_ids[first_read_idx]); 

        if (read_obj->hasConstrainedCluster()) {
            con_c_id = read_obj->getConstrainedClusterID();
            
            // Only enforce if the PARTNER is also in this column!
            if (cluster_id_to_read_index_map.count(con_c_id) > 0) {
                dependent = std::max(c_id, con_c_id);
                anchor = std::min(c_id, con_c_id);
                read_cluster_constraints[dependent] = anchor;
            }
        }
    }

	/**
	 * Figuring out the positions which should be subject to Gray Code ordering.
	 * We also find the constrained positions mapping for updating the binary vector later.
	 * 
	 * Example:
	 * Let the read clusters be {1,2,5,6,8,9,14} and let the constraints be 6->2 and 9->8.
	 * Then the free positions will become {0,1,2,4,6} 
	 * 
	 * We perform Gray code ordering for {0, 1, 2, 4, 6}.
	 * 
	 * And the constrained position mapping will be:
	 *   1 -> 3 => index 1 in read_clusters (which is cluster ID 2) maps to index 3 (which is cluster ID 6)
	 *   4 -> 5 => index 4 in read_clusters (which is cluster ID 8) maps to index 5 (which is cluster ID 9)
	 */
	uint32_t pos;
	uint32_t pos_c_id;
	uint32_t j;
	std::vector<uint32_t> free_read_cluster_positions;		// the position of the free read clusters in read_cluster_ids
	std::vector<uint32_t> num_reads_per_free_read_cluster_positions;	// tracks the number of reads associated with the read cluster in free_read_cluster_positions;
	bool found;
	for (uint32_t i = 0; i < read_cluster_ids.size(); i++) {
		c_id = read_cluster_ids[i];
		if (read_cluster_constraints.count(c_id) == 0) {
			/**
			 * This cluster ID is either untagged (and hence unconstrained)
			 * Or this is the min of the constrained pair (given as max -> min mapping). 
			 * Hence only the min will be varied in Gray Code.
			 */
			free_read_cluster_positions.push_back(i);
			if (cluster_id_to_read_index_map.count(c_id) == 0) { num_reads_per_free_read_cluster_positions.push_back(1); } // cluster is actually just a single read.
			else { num_reads_per_free_read_cluster_positions.push_back(cluster_id_to_read_index_map.count(c_id)); }
		}
		else {
			con_c_id = read_cluster_constraints.at(c_id); // con_c_id is the cluster ID which is included in free_read_cluster_positions
			/**
			 * Now we find the position of con_c_id in the free_read_cluster_positions vector.
			 * We know that con_c_id must be already included in free_read_cluster_positions since it is the min of the constrained pair.
			 */
			found = false;
			for (j = 0; j < free_read_cluster_positions.size(); j++) {
				pos = free_read_cluster_positions[j];
				pos_c_id = read_cluster_ids[pos];
				if (pos_c_id == con_c_id) {
					// Found the position of representative cluster ID in free_read_cluster_positions vector.
					// Now we can store the mapping.
					constrained_position_map[pos] = i;
					// We will also update the num_reads_per_free_read_cluster_positions
					num_reads_per_free_read_cluster_positions[j] += cluster_id_to_read_index_map.count(c_id);
					found = true;
					break;
				}
			}
			if (!found) {
				/**
				 * The rep_c_id was not found in free_read_cluster_positions.
				 * This indicates that the position was heterozygous or an SV which was not considered for phasing.
				 * Hence we can consider this cluster as free to vary in Gray Code ordering.
				 * 
				 */
				assert(false);
				free_read_cluster_positions.push_back(i);
				// constraine_position_map[pos] does not exist.
			}
		}
	}

	/**
	 * Sorting the free_read_cluster_positions based on how many reads are associated with the cluster given in num_reads_per_free_read_cluster_positions.
	 * If the cluster is constrained, then the reads in the constrained cluster are also taken into consideration.
	 */
    std::vector<std::pair<uint32_t, uint32_t>> sort_helper;
    for(uint32_t k = 0; k < free_read_cluster_positions.size(); k++) {
        sort_helper.push_back({num_reads_per_free_read_cluster_positions[k], free_read_cluster_positions[k]});
    }
    std::sort(sort_helper.begin(), sort_helper.end());
    for(const auto& p : sort_helper) {
        sorted_free_read_cluster_positions.push_back(p.second);
    }
	/**
	 * Finding the clusters from the next variant position
	 */
	for (const auto& next_read_id : next_read_ids) {
        read_obj = set->get(next_read_id);
        c_id = (!read_obj->isClustered()) ? next_read_id : read_obj->getClusterID();
        next_read_cluster_ids.push_back(c_id);
    }
	/**
	 * The bipartition logic of backward compatibility requires both
	 * read_cluster_ids and next_read_cluster_ids to be sorted.
	 */
	std::sort(next_read_cluster_ids.begin(), next_read_cluster_ids.end());
    // Remove duplicates from next_read_cluster_ids after sorting
    auto last = std::unique(next_read_cluster_ids.begin(), next_read_cluster_ids.end());
    next_read_cluster_ids.erase(last, next_read_cluster_ids.end());

	precompute_bipartition(next_read_cluster_ids);
}

unique_ptr<BipartitionIterator> Column::get_iterator(ReadSet* set) {
	return unique_ptr<BipartitionIterator>(new BipartitionIterator(this, set));
}


const vector<uint32_t> * Column::get_read_cluster_ids() const {
	return &(this->read_cluster_ids);
}

const vector<uint32_t> * Column::get_read_ids() const {
	return &(this->read_ids);
}

const unordered_map<uint32_t, vector<uint32_t>> * Column::get_cluster_id_to_read_index_map() const {
	return &(this->cluster_id_to_read_index_map);
}


const vector<uint32_t> * Column::get_sorted_free_read_cluster_positions() const {
	return &(this->sorted_free_read_cluster_positions);
}

const unordered_map<uint32_t, uint32_t> * Column::get_constrained_position_map() const {
	return &(this->constrained_position_map);
}


/**
 * The number of bipartitions is 2^{number of read clusters - number of constrained cluster pairs}
 * We check for the number of constrained pairs using read_cluster_contraints.
 * 
 * @see BipartitionIterator's constructor (at the end)
 */
uint32_t Column::get_num_bipartition() const {
	uint32_t num_free_read_cluster = read_cluster_ids.size() - constrained_position_map.size();

	// Safety check for shift overflow (if num_free_read_cluster_positions >= 32)
    if (num_free_read_cluster >= 32) {
         throw std::overflow_error("Too many partitions for uint32_t");
	}
	
	return (1 << num_free_read_cluster);
}

/**
 * Precomputes the compatible bipartitions in the left column (at position index)
 * Since we already have the read clusters, we know which clusters are common between the two columns.
 * The compatible bipartitions simply fix the bits for the common clusters as per the b_index of the right column (at position index + 1).
 * 
 * So in this function, we do a precomputation of all the compatible bipartitions based on the read clusters that are unique to the left column (at position index).
 * Then during runtime, we just need to add the fixed bits from b_index of the right column to these precomputed bipartitions.
 */
void Column::precompute_bipartition(std::vector<uint32_t>& next_read_cluster_ids) {
	
	// initialize cached_bipartitions
	cached_bipartitions = {0};

	/**
	 * Create a map between Cluster ID and its index in sorted_free_read_cluster_positions.
	 * Since the position of the cluster determines where it is in the Gray Code ordering.
	 * This will be used to define the masks for next_read_cluster_masks.
	 */
	std::unordered_map<uint32_t, uint32_t> cluster_id_to_graycode_index_map;
	uint32_t count = 0;
	uint32_t cluster_id;
	for (uint32_t cluster_pos: sorted_free_read_cluster_positions) {
		cluster_id = read_cluster_ids[cluster_pos];
		cluster_id_to_graycode_index_map[cluster_id] = count;
		count++;
	}

	/**
	 * initialize masks for all the read clusters in the next column.
	 * some of them will remain 0 if they are unique to the next column.
	 */
	next_read_cluster_masks.resize(next_read_cluster_ids.size(), 0);
	
	/**
	 * Example to explain the logic:
	 * Let read_ids = {1,2,3,4,5,6,7,8,9,10} and next_read_ids = {2,3,5,6,8}
	 * Let the read clusters be 3, 4 as 3 and 7, 8 as 7.
	 * So read_cluster_ids = {1,2,3,5,6,7,9,10} and next_read_cluster_ids = {2,3,6,7}
	 * We want to precompute the compatible bipartitions in the left column since we know that clusters 1,5,9,10 are unique to the left column.
	 * So we iterate through next_read_cluster_ids and read_cluster_ids simultaneously.
	 */
	count = 0;	// count keeps track of the position in read_cluster_ids
	uint32_t next_id;
	uint32_t current_size;
	uint32_t bit_val;
	for (uint32_t i = 0; i < next_read_cluster_ids.size(); i++) {
		next_id = next_read_cluster_ids[i];

		/**
		 * The while loop iterates through all the read clusters in the left column (at position index)
		 * that are less than the current read cluster in the right column (at position index + 1).
		 * For these read clusters, we can have both possibilities (0 and 1) since they are not present
		 * in the right column.
		 * 
		 * Hence we double the size of compatible_bipartition and add the corresponding values.
		 * 
		 * From the example above, the whie loop processing the clusters 1 and 5.
		 */
		if (count < read_cluster_ids.size()) {
            while ((count < read_cluster_ids.size()) && (next_id > read_cluster_ids[count])) {
				if (cluster_id_to_graycode_index_map.count(read_cluster_ids[count]) == 0) {
					// this cluster is not part of the Gray Code
					count++;
					continue;
				}
                current_size = cached_bipartitions.size();
                cached_bipartitions.resize(2 * current_size);
                
                // Pre-calculate the bit value
                bit_val = (1 << cluster_id_to_graycode_index_map[read_cluster_ids[count]]); 

                for (uint32_t j = 0; j < current_size; j++) {
                    cached_bipartitions[current_size + j] = cached_bipartitions[j] + bit_val;
                }
				count++;
            }
        }
		/**
		 * When we exhaust all the clusters less than ri (in the while loop above),
		 * we check if the ri exists in the read_cluster_ids (left column).
		 * 
		 * It has to! Either ri is present in both columns or 
		 * it is a new cluster that was found in the next position.
		 * 
		 * If it is the later, then count == read_cluster_ids.size() and we exit the for loop.
		 * 
		 * Here I am storing the position of the shared cluster in the left column.
		 * This will be useful during runtime to quickly determine the compatible bipartitions.
		 * 
		 * From the example above, the for loop processing the clusters 2,3,6,7.
		 */
		if (count < read_cluster_ids.size() && next_id == read_cluster_ids.at(count)) {
			// Storing the position of the shared cluster in the left column
			if (cluster_id_to_graycode_index_map.count(read_cluster_ids[count]) != 0) {
				// checking if this cluster is actually part of gray code.
				next_read_cluster_masks[i] = (1 << cluster_id_to_graycode_index_map[read_cluster_ids[count]]);
			}
			count++;
		}
	}

	/**
	 * Here we process all the new clusters are exclusively in the left column (at position index).
	 * These are the clusters that are present in read_cluster_ids but not in next_read_cluster_ids even after the above for loop.
	 * 
	 * From the example above, the while loop processing the clusters 9,10.
	 */
	for (uint32_t i = count; i < read_cluster_ids.size(); i++) {
		if (cluster_id_to_graycode_index_map.count(read_cluster_ids[i]) == 0) {
			// this cluster is not part of the Gray Code
			continue;
		}
        current_size = cached_bipartitions.size();
        cached_bipartitions.resize(2 * current_size);
        
        bit_val = (1 << cluster_id_to_graycode_index_map[read_cluster_ids[i]]);

        for (uint32_t j = 0; j < current_size; j++) {
            cached_bipartitions[current_size + j] = cached_bipartitions[j] + bit_val;
        }
    }
}

/**
 * Using the precomputed bipartition, find the compatible bipartitions in the left column (at position index)
 * for the given read_cluster_bit_representation of the right column (at position index + 1).
 * 
 * NOTE: The read_cluster_bit_representation contains the bit information of which bipartition each read cluster is in.
 */
void Column::get_backward_compatible_bipartitions(uint32_t read_cluster_bit_representation, vector<uint32_t>& result) const {
	const uint32_t n = cached_bipartitions.size();
    if (result.size() != n) result.resize(n);	
	
	/**
	 * base keeps track of the value to be added to all compatible bipartitions.
	 * This will be based on the read clusters that are common between the two columns.
	 * In the example below, read clusters 2,3,6,7 are common between the two columns.
	 * So their bit assignments should be the same in both columns.
	 */
	uint32_t base = 0;

	for (uint32_t mask : next_read_cluster_masks) {
        /**
		 * if the current cluster has a bit value of 1 in b_index of the right column,
		 * then we need to set the corresponding bit in base.
		 * If the cluster is not present in the right column, then its mask will be 0
		 * and we skip it.
		 */
        if (read_cluster_bit_representation & 1) {
            base += mask;
        }
        
        // Shift b_index to the next bit
        read_cluster_bit_representation >>= 1;
    }

    // Applying the base to all compatible bipartitions
	if (base > 0) {
        for (size_t i = 0; i < result.size(); i++) {
            result[i] += cached_bipartitions[i] + base;
        }
    }
}
