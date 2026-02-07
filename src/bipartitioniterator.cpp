// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#include <cassert>
#include "column.h"
#include "bipartitioniterator.h"
#include "readset.h"
#include "math.h"

BipartitionIterator::BipartitionIterator(Column* parent, ReadSet* set) {
	assert(parent != 0);
	this->parent = parent;

	/**
	 * Figuring out the positions which should be subject to Gray Code ordering.
	 * We also find the constrained positions mapping for updating the binary vector later.
	 * 
	 * Example:
	 * Let the read clusters be {1,2,5,6,8,9,14} and let the constraints be 6->2 and 9->8.
	 * Then the free positions will become {0,1,2,4,6} 
	 * 
	 * We perform Gray code ordering for {0, 1, 2, 4} and fix 6 to be in biparition 0
	 * since bipartition is complementary (discussed below right above the new GrayCodes statement)
	 * 
	 * And the constrained position mapping will be:
	 *   1 -> 3 => index 1 in read_clusters (which is cluster ID 2) maps to index 3 (which is cluster ID 6)
	 *   4 -> 5 => index 4 in read_clusters (which is cluster ID 8) maps to index 5 (which is cluster ID 9)
	 */
	const std::unordered_map<uint32_t, uint32_t> read_cluster_constraints = *(parent->get_read_cluster_constraints_map());
	for (uint32_t i = 0; i < parent->get_read_cluster_ids()->size(); i++) {
		uint32_t c_id = parent->get_read_cluster_ids()->at(i);
		if (read_cluster_constraints.count(c_id) == 0) {
			/**
			 * This cluster ID is either untagged (and hence unconstrained)
			 * Or this is the min of the constrained pair (given as max -> min mapping). 
			 * Hence only the min will be varied in Gray Code.
			 */
			free_positions.push_back(i);
		}
		else {
			uint32_t rep_c_id = read_cluster_constraints.at(c_id); // rep_c_id is the cluster ID which is included in free_positions
			/**
			 * Now we find the position of rep_c_id in the free_positions vector.
			 * We know that rep_c_id must be already included in free_positions since it is the min of the constrained pair.
			 */
			bool found = false;
			for (uint32_t j = 0; j < free_positions.size(); j++) {
				uint32_t pos = free_positions[j];
				uint32_t pos_c_id = parent->get_read_cluster_ids()->at(pos);
				if (pos_c_id == rep_c_id) {
					// Found the position of representative cluster ID in free_positions vector.
					// Now we can store the mapping.
					constrained_position_map[pos] = i;
					found = true;
					break;
				}
			}
			if (!found) {
				/**
				 * The rep_c_id was not found in free_positions.
				 * This indicates that the position was heterozygous or an SV which was not considered for phasing.
				 * Hence we can consider this cluster as free to vary in Gray Code ordering.
				 */
				free_positions.push_back(i);
				// constraine_position_map[pos] does not exist.
			}
		}
	}
	/**
	 * Initialize the bit representation to 0 indicating all clusters
	 * are in bipartition 0.
	 * 
	 * The non-represented cluster IDs from constraint pairs will 
	 * then be flipped since the pair can not be in the same bipartition.
	 * 
	 * So for the read clusters {1,2,5,6,8,9,14} and constraints 6->2 and 9->8 (as given above in the example),
	 * will have the following initial representation:
	 * 
	 * clusters -> 14 | 9 | 8 | 6 | 5 | 2 | 1
	 * bits -> .... 0 | 1 | 0 | 1 | 0 | 0 | 0  -> 40
	 */ 
	this->read_cluster_bit_representation = 0;
	for (auto pos : constrained_position_map) {
		// Setting the non-representative constrained clusters to true initially.
		// since these are the only positions which are initially true, we can set b_index accordingly.
		int mask = 1 << pos.second;
		this->read_cluster_bit_representation |= mask;
	}

	/**
	 * We Gray code order with free_position.size() - 1 because 
	 * the read bipartition assignment in complementary.
	 * 
	 * Read 1, 2 in bipartition 1 and Read 3, 4 in bipartition 2
	 * is equivalent to
	 * Read 1, 2 in bipartition 2 and Read 3, 4 in bipartition 1
	 * 
	 * So, with this implementation, from the above example, 
	 * we vary clusters 1, 2, 5, and 8. Cluster 14 is fixed to 0.
	 */
	this->graycodes = new GrayCodes(free_positions.size() - 1);
}


BipartitionIterator::~BipartitionIterator() {
	delete graycodes;
}


bool BipartitionIterator::has_next() const{
	return graycodes->has_next();
}

/**
 * Advance to the next bipartition index using Gray Code ordering.
 * Since we update only the free positions in Gray Code ordering, we need to translate the changed bit to which clusters are changed.
 * The following things need to happen:
 * 1. The bit flipped in the Gray Code (given by free_positions) needs to translated to which cluster (given by binary_vector) is changed.
 * 2. If the cluster changed is constrained with another cluster, then that other cluster also needs to be updated in binary_vector.
 * 3. The b_index (which is based on binary_vector) needs to be updated accordingly.
 */
void BipartitionIterator::advance(int* cluster_bit_changed) {
	
	assert(graycodes->has_next());
	/**
	 * graycode_bit_changed tracks which bit in the Gray Code ordering was flipped.
	 * Consisdering the example from the Constructor, if the free_positions = {1,2,5,8,14},
	 * then graycode_bit_changed tells which bit was flipped was last iteration.
	 * 
	 * graycode_bit_changed = x => xth index of free_positions was flipped (0->1, 1->2, 2->5, 3->8, 4->14)
	 * 
	 * graycode_bit_changed = -1 => initialisation step and all the bits are 0.
	 */
	int graycode_bit_changed = -1;
	/**
	 * bipartition_index stores the current state of the Gray Code ordering (after doing the flipping).
	 * Let's say the state (before executing get_next()) was 00000.
	 * After executing get_next(), the new ordering becomes 00001.
	 *   - graycode_bit_changed = 0 (since the first bit was flipped)
	 *   - graycode_binaryindex = 1 (which is the numerical representation of 00001)
	 */
	this->bipartition_index = graycodes->get_next(&graycode_bit_changed);
	
	// finding which bit in binary_vector is changed based on the graycode_bit_changed.
	if (cluster_bit_changed != 0) {
		if (graycode_bit_changed == -1) {
			*cluster_bit_changed = -1;
		}
		else {
			// bit_changed gets the value of the cluster id whose bit was flipped.
			*cluster_bit_changed = this->free_positions[graycode_bit_changed]; // this is the position in binary_vector which is changed.
		}
	}

	// now we check if graycode_bit_changed corresponds to a constrained cluster.
	if (graycode_bit_changed != -1 && constrained_position_map.count(*cluster_bit_changed) > 0) {
		// The changed cluster is constrained with another cluster.
		// constrained_position gives us the positions of the non-representative cluster in binary_vector.
		uint32_t constrained_position = constrained_position_map.at(*cluster_bit_changed);
		// Updating b_index using masks.
		int mask1 = 1 << *cluster_bit_changed; // mask for the representative cluster.
		this->read_cluster_bit_representation ^= mask1; // XOR operation to flip the bit.
		int mask2 = 1 << constrained_position; // mask for the non-representative cluster.
		this->read_cluster_bit_representation ^= mask2; // XOR operation to flip the bit.
	} else if (graycode_bit_changed != -1) {
		// Updating b_index using masks.
		int mask = 1 << *cluster_bit_changed;
		this->read_cluster_bit_representation ^= mask;
	} else {
		// Initialisation step where all the bits are set according to graycode_binaryvector.
		// since b_index is already set correctly in the constructor, we do not need to update it here.
	}
}

Column* BipartitionIterator::get_parent_column() const {
	return parent;
}

uint32_t BipartitionIterator::get_bipartition_index() const {
	assert(this->bipartition_index < num_bipartitions());
	return this->bipartition_index;
}

uint32_t BipartitionIterator::get_read_cluster_bit_representation() const {
	return this->read_cluster_bit_representation;
}

uint32_t BipartitionIterator::num_bipartitions() const {
	if (free_positions.empty()) {
		throw std::runtime_error("No read clusters found! Terminating!");
	}
    
    // Check if size is 1 (only 1 cluster -> 1 partition: 0)
    if (free_positions.size() == 1) return 1;

    uint32_t num_free_positions = free_positions.size() - 1;
    
    // Safety check for shift overflow (if num_free_positions >= 32)
    if (num_free_positions >= 32) {
         throw std::overflow_error("Too many partitions for uint32_t");
    }

    return (1 << num_free_positions);
}

bool BipartitionIterator::is_clustered_bit(uint32_t cluster_bit_changed) const {
	return parent->get_cluster_id_to_read_index_map()->count(cluster_bit_changed) > 0;
}

bool BipartitionIterator::has_constrained_bit(uint32_t cluster_bit_changed) const {
	return constrained_position_map.count(cluster_bit_changed) > 0;
}

uint32_t BipartitionIterator::get_constrained_bit(uint32_t cluster_bit_changed) const {
	assert(is_clustered_bit(cluster_bit_changed));
	assert(has_constrained_bit(cluster_bit_changed));
	return constrained_position_map.at(cluster_bit_changed);
}

std::vector<uint32_t> BipartitionIterator::get_read_index_from_cluster_id(uint32_t cluster_id) const {
	const std::unordered_map<uint32_t, std::vector<uint32_t>> cluster_id_to_read_index_map = *(parent->get_cluster_id_to_read_index_map());
	return cluster_id_to_read_index_map.at(cluster_id);
}

std::unordered_map<uint32_t, bool> BipartitionIterator::get_changed_reads(uint32_t cluster_bit_changed) const {
	std::unordered_map<uint32_t, bool> changed_bits;
	bool new_bit = (read_cluster_bit_representation >> cluster_bit_changed) & 1;
	uint32_t flipped_cluster_id = parent->get_read_cluster_ids()->at(cluster_bit_changed);
	if (!is_clustered_bit(cluster_bit_changed)) {
		// The bit did not correspond to a cluster.
		changed_bits[flipped_cluster_id] = new_bit; // cluster ID is the read ID.
		return changed_bits;
	}
	/**
	 * The bit changed corresponds to a cluster.
	 * Adding the information of the read_ids corresponding to this cluster.
	 */
	std::vector<uint32_t> read_indices_from_clusters = get_read_index_from_cluster_id(flipped_cluster_id);
	for (auto read_index: read_indices_from_clusters) {
		changed_bits[read_index] = new_bit;
	}
	/**
	 * If there is a constrained cluster, add the reads id from that cluster.
	 */
	if (has_constrained_bit(cluster_bit_changed)) {
		uint32_t constrained_flipped_cluster_id = parent->get_read_cluster_ids()->at(constrained_position_map.at(cluster_bit_changed));
		std::vector<uint32_t> constrained_read_index_from_clusters = get_read_index_from_cluster_id(constrained_flipped_cluster_id);
		for (auto read_index: constrained_read_index_from_clusters) {
			changed_bits[read_index] = !new_bit;
		}
	}
	return changed_bits;
}