// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#include <cassert>
#include "column.h"
#include "columnindexingiterator.h"
#include "readset.h"
#include "math.h"

ColumnIndexingIterator::ColumnIndexingIterator(Column* parent, ReadSet* set) {
	assert(parent != 0);
	this->parent = parent;

	/**
	 * Figuring out the positions which should be subject to Gray Code ordering.
	 * We also find the constrained positions mapping for updating the binary vector later.
	 * 
	 * Example:
	 * Let the read clusters be {1,2,5,6,8,9,14} and let the constraints be 6->2 and 9->8.
	 * Then the free positions will become {0,1,2,4,6} (these are the positions which will be varied in Gray Code)
	 * And the constrained position mapping will be:
	 *   1 -> 3 => index 1 in binary_vector (which is cluster ID 2) maps to index 3 in binary_vector
	 *   4 -> 5 => index 4 in binary_vector (which is cluster ID 8) maps to index 5 in binary_vector
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
			uint32_t rep_c_id = read_cluster_constraints.at(c_id); // this is the cluster ID which is included in free_positions
			/**
			 * Now we find the position of rep_c_id in the free_positions vector.
			 * We know that rep_c_id must be already included in free_positions since it is the min of the constrained pair.
			 */
			for (uint32_t j = 0; j < free_positions.size(); j++) {
				uint32_t pos = free_positions[j];
				uint32_t pos_c_id = parent->get_read_cluster_ids()->at(pos);
				if (pos_c_id == rep_c_id) {
					// Found the position of representative cluster ID in free_positions vector.
					// Now we can store the mapping.
					constrained_position_map[pos] = i;
					break;
				}
			}
		}
	}
	/**
	 * Initialize the binary vector to false for all clusters represented in free_positions.
	 * The non-represented cluster IDs will get true as their initial bipartition assignment.
	 * So for the read clusters {1,2,5,6,8,9,14} and constraints 6->2 and 9->8 (as given above in the example),
	 * we get the initial binary vector as {false, false, false, true, false, true, false}
	 */ 
	this->b_index = 0;
	for (auto pos : constrained_position_map) {
		// Setting the non-representative constrained clusters to true initially.
		// since these are the only positions which are initially true, we can set b_index accordingly.
		int mask = 1 << pos.second;
		this->b_index |= mask;
	}


	// The Gray Code ordering is now done only on the positions whose bipartitions are unknown
	this->graycodes = new GrayCodes(free_positions.size());
}


ColumnIndexingIterator::~ColumnIndexingIterator() {
	delete graycodes;
}


bool ColumnIndexingIterator::has_next() {
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
void ColumnIndexingIterator::advance(int* bit_changed) {
	
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
	 * graycode_binaryindex gives the current state of the Gray Code ordering (after doing the flipping).
	 * Let's say the state (before executing get_next()) was 00000.
	 * After executing get_next(), the new ordering becomes 00001.
	 *   - graycode_bit_changed = 0 (since the first bit was flipped)
	 *   - graycode_binaryindex = 1 (which is the numerical representation of 00001)
	 */
	uint32_t graycode_binaryindex = graycodes->get_next(&graycode_bit_changed);
	
	// finding which bit in binary_vector is changed based on the graycode_bit_changed.
	if (bit_changed != 0) {
		if (graycode_bit_changed == -1) {
			*bit_changed = -1;
		}
		else {
			// bit_changed gets the value of the cluster id whose bit was flipped.
			*bit_changed = this->free_positions[graycode_bit_changed]; // this is the position in binary_vector which is changed.
		}
	}

	// now we check if graycode_bit_changed corresponds to a constrained cluster.
	if (graycode_bit_changed != -1 && constrained_position_map.count(*bit_changed) > 0) {
		// The changed cluster is constrained with another cluster.
		// constrained_position gives us the positions of the non-representative cluster in binary_vector.
		uint32_t constrained_position = constrained_position_map.at(*bit_changed);
		// Updating b_index using masks.
		int mask1 = 1 << *bit_changed; // mask for the representative cluster.
		this->b_index ^= mask1; // XOR operation to flip the bit.
		int mask2 = 1 << constrained_position; // mask for the non-representative cluster.
		this->b_index ^= mask2; // XOR operation to flip the bit.
	} else if (graycode_bit_changed != -1) {
		// Updating b_index using masks.
		int mask = 1 << *bit_changed;
		this->b_index ^= mask;
	} else {
		// Initialisation step where all the bits are set according to graycode_binaryvector.
		// since b_index is already set correctly in the constructor, we do not need to update it here.
	}
}

uint32_t ColumnIndexingIterator::get_b_index() {
	return this->b_index;
}

bool ColumnIndexingIterator::is_clustered_bit(uint32_t bit_changed) {
	return parent->get_cluster_id_to_read_ids_map()->count(parent->get_read_cluster_ids()->at(bit_changed)) > 0;
}

bool ColumnIndexingIterator::has_constrained_bit(uint32_t bit_changed) {
	return constrained_position_map.count(bit_changed) > 0;
}

uint32_t ColumnIndexingIterator::get_constrained_bit(uint32_t bit_changed) {
	assert(is_clustered_bit(bit_changed));
	return constrained_position_map.at(bit_changed);
}

std::vector<uint32_t> ColumnIndexingIterator::get_reads_from_cluster_id(uint32_t cluster_id) {
	const std::unordered_map<uint32_t, std::vector<uint32_t>> cluster_id_to_read_ids_map = *(parent->get_cluster_id_to_read_ids_map());
	return cluster_id_to_read_ids_map.at(cluster_id);
}

std::unordered_map<uint32_t, bool> ColumnIndexingIterator::get_changed_bits(uint32_t bit_changed) {
	std::unordered_map<uint32_t, bool> changed_bits;
	bool new_bit = (b_index >> bit_changed) & 1;
	uint32_t flipped_cluster_id = parent->get_read_cluster_ids()->at(bit_changed);
	if (!is_clustered_bit(bit_changed)) {
		// The bit did not correspond to a cluster.
		changed_bits[flipped_cluster_id] = new_bit; // cluster ID is the read ID.
		return changed_bits;
	}
	/**
	 * The bit changed corresponds to a cluster.
	 * Adding the information of the read_ids corresponding to this cluster.
	 */
	std::vector<uint32_t> read_ids_from_clusters = get_reads_from_cluster_id(flipped_cluster_id);
	for (auto read_id: read_ids_from_clusters) {
		changed_bits[read_id] = new_bit;
	}
	/**
	 * If there is a constrained cluster, add the reads id from that cluster.
	 */
	if (has_constrained_bit(bit_changed)) {
		uint32_t constrained_flipped_cluster_id = parent->get_read_cluster_ids()->at(constrained_position_map[bit_changed]);
		std::vector<uint32_t> constrained_read_ids_from_clusters = get_reads_from_cluster_id(constrained_flipped_cluster_id);
		for (auto read_id: constrained_read_ids_from_clusters) {
			changed_bits[read_id] = !new_bit;
		}
	}
	return changed_bits;
}