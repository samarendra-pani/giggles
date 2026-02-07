// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#include <cassert>
#include <math.h>
#include <algorithm>

#include "column.h"

using namespace std;

Column::Column(const std::vector<uint32_t>& read_ids, const std::vector<uint32_t>& next_read_ids, ReadSet* set): 
	read_ids(read_ids) {
	
	std::vector<uint32_t> next_read_cluster_ids;
	/**
	 * Finding the read clusters from the phasing (and possibly other clustering later)
	 */
	uint32_t count = 0;
	for (const auto& current_read_id : read_ids) {
        Read* read_obj = set->get(current_read_id);

        if (!read_obj->isClustered()) {
			// Read is not clustered. So we use it as its own cluster.
            read_cluster_ids.push_back(current_read_id);
        } else {
            uint32_t cluster_id = read_obj->getClusterID();
            
            // Only add unique cluster IDs to the list
            if (cluster_id_to_read_index_map.count(cluster_id) == 0) {
                read_cluster_ids.push_back(cluster_id);
            }
            cluster_id_to_read_index_map[cluster_id].push_back(count);
        }
		count ++;
    }
	/**
	 * Finding the read cluster constraints from the cluster id in the cluster_id_to_read_index_map
	 */
	for (const auto& pair1 : cluster_id_to_read_index_map) {
		uint32_t c_id1 = pair1.first;
		Read* read_obj1 = set->get(c_id1);
		if (read_obj1->hasConstrainedCluster()) {
			uint32_t c_id2 = read_obj1->getConstrainedClusterID();
			// FORCE DIRECTION: Max -> Min
			uint32_t dependent = std::max(c_id1, c_id2);
			uint32_t anchor = std::min(c_id1, c_id2);

			// This automatically handles the "double counting" check since we always store max -> min
			read_cluster_constraints[dependent] = anchor;
		}
	}
	/**
	 * Finding the clusters from the next variant position
	 */
	for (const auto& next_read_id : next_read_ids) {
        Read* read_obj = set->get(next_read_id);
        uint32_t id_to_add = (!read_obj->isClustered()) ? next_read_id : read_obj->getClusterID();
        next_read_cluster_ids.push_back(id_to_add);
    }
	/**
	 * The bipartition logic of backward compatibility requires both
	 * read_cluster_ids and next_read_cluster_ids to be sorted.
	 */
	std::sort(read_cluster_ids.begin(), read_cluster_ids.end());
    std::sort(next_read_cluster_ids.begin(), next_read_cluster_ids.end());
    // Remove duplicates from next_read_cluster_ids after sorting
    auto last = std::unique(next_read_cluster_ids.begin(), next_read_cluster_ids.end());
    next_read_cluster_ids.erase(last, next_read_cluster_ids.end());

	precompute_bipartition(next_read_cluster_ids);
}

unique_ptr<BipartitionIterator> Column::get_iterator(ReadSet* set) {
	return unique_ptr<BipartitionIterator>(new BipartitionIterator(this, set));
}


vector<uint32_t> * Column::get_read_cluster_ids() {
	return &(this->read_cluster_ids);
}

vector<uint32_t> * Column::get_read_ids() {
	return &(this->read_ids);
}

unordered_map<uint32_t, vector<uint32_t>> * Column::get_cluster_id_to_read_index_map() {
	return &(this->cluster_id_to_read_index_map);
}


unordered_map<uint32_t, uint32_t> * Column::get_read_cluster_constraints_map() {
	return &(this->read_cluster_constraints);
}
/**
 * Precomputes the compatible bipartitions in the left column (at position index)
 * Since we already have the read clusters, we know which clusters are common between the two columns.
 * The compatible bipartitions simply fix the bits for the common clusters as per the b_index of the right column (at position index + 1).
 * 
 * So in this function, we do a precomputation of all the compatible bipartitions based on the read clusters that are unique to the left column (at position index).
 * Then during runtime, we just need to add the fixed bits from b_index of the right column to these precomputed bipartitions.
 */
void Column::precompute_bipartition(std::vector<uint32_t> next_read_cluster_ids) {
	
	// initialize cached_bipartitions
	cached_bipartitions = {0};

	/**
	 * initialize masks for all the read clusters in the next column.
	 * some of them will remain 0 if they are unique to the next column.
	 */
	next_read_cluster_masks.resize(next_read_cluster_ids.size(), 0);
	
	// count keeps track of the position in read_cluster_ids
    uint32_t count = 0;
	
	/**
	 * Requirement: read_cluster_ids and next_read_cluster_ids should be sorted in ascending order.
	 */
	/**
	 * Example to explain the logic:
	 * Let read_ids = {1,2,3,4,5,6,7,8,9,10} and next_read_ids = {2,3,5,6,8}
	 * Let the read clusters be 3, 4 as 3 and 7, 8 as 7.
	 * So read_cluster_ids = {1,2,3,5,6,7,9,10} and next_read_cluster_ids = {2,3,6,7}
	 * We want to precompute the compatible bipartitions in the left column since we know that clusters 1,5,9,10 are unique to the left column.
	 * So we iterate through next_read_cluster_ids and read_cluster_ids simultaneously.
	 */
	for (uint32_t i = 0; i < next_read_cluster_ids.size(); i++) {
		uint32_t next_id = next_read_cluster_ids.at(i);

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
            while ((count < read_cluster_ids.size()) && (next_id > read_cluster_ids.at(count))) {
                uint32_t current_size = cached_bipartitions.size();
                cached_bipartitions.resize(2 * current_size);
                
                // Pre-calculate the bit value
                uint32_t bit_val = (1 << count); 

                for (uint32_t j = 0; j < current_size; j++) {
                    cached_bipartitions.at(current_size + j) = cached_bipartitions.at(j) + bit_val;
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
			next_read_cluster_masks.at(i) = (1 << count);
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
        uint32_t current_size = cached_bipartitions.size();
        cached_bipartitions.resize(2 * current_size);
        
        uint32_t bit_val = (1 << i);

        for (uint32_t j = 0; j < current_size; j++) {
            cached_bipartitions.at(current_size + j) = cached_bipartitions.at(j) + bit_val;
        }
    }
}

/*
 * Using the precomputed bipartition, find the compatible bipartitions in the left column (at position index)
 * for the given b_index of the right column (at position index + 1).
 */
vector<uint32_t> Column::get_backward_compatible_bipartitions(uint32_t b_index) {
	vector<uint32_t> result = cached_bipartitions;	
	
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
        if (b_index & 1) {
            base += mask;
        }
        
        // Shift b_index to the next bit
        b_index >>= 1;
    }

    // Applying the base to all compatible bipartitions
	if (base > 0) {
        for (size_t i = 0; i < result.size(); i++) {
            result[i] += base;
        }
    }

	return result;
}
