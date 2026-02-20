#include <cassert>

#include "bipartitioniterator.h"

BipartitionIterator::BipartitionIterator(Column* parent, ReadSet* set) {
	assert(parent != 0);
	this->parent = parent;

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
	this->bipartition_index = 0;
	uint32_t mask;
	for (auto pos : *(parent->get_constrained_position_map())) {
		// Setting the non-representative constrained clusters to true initially.
		// since these are the only positions which are initially true, we can set b_index accordingly.
		mask = 1 << pos.second;
		this->read_cluster_bit_representation |= mask;
	}
	this->graycodes = new GrayCodes(parent->get_sorted_free_read_cluster_positions()->size());
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
 * 1. The bit flipped in the Gray Code (given by free_read_cluster_positions) needs to translated to which cluster (given by binary_vector) is changed.
 * 2. If the cluster changed is constrained with another cluster, then that other cluster also needs to be updated in binary_vector.
 * 3. The b_index (which is based on binary_vector) needs to be updated accordingly.
 */
void BipartitionIterator::advance(int* cluster_bit_changed) {
	
	assert(graycodes->has_next());
	/**
	 * graycode_bit_changed tracks which bit in the Gray Code ordering was flipped.
	 * Consisdering the example from the Constructor, if the free_read_cluster_positions = {1,2,5,8,14},
	 * then graycode_bit_changed tells which bit was flipped was last iteration.
	 * 
	 * graycode_bit_changed = x => xth index of free_read_cluster_positions was flipped (0->1, 1->2, 2->5, 3->8, 4->14)
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
			// this is the bit position in read_cluster_bit_representation which is changed.
			*cluster_bit_changed = this->parent->get_sorted_free_read_cluster_positions()->at(graycode_bit_changed);
		}
	}

	// now we check if graycode_bit_changed corresponds to a constrained cluster.
	if (graycode_bit_changed != -1 && this->parent->get_constrained_position_map()->count(*cluster_bit_changed) > 0) {
		// The changed cluster is constrained with another cluster.
		// constrained_position gives us the positions of the non-representative cluster in binary_vector.
		uint32_t constrained_position = this->parent->get_constrained_position_map()->at(*cluster_bit_changed);
		// Updating b_index using masks.
		uint32_t mask1 = 1 << *cluster_bit_changed; // mask for the representative cluster.
		this->read_cluster_bit_representation ^= mask1; // XOR operation to flip the bit.
		uint32_t mask2 = 1 << constrained_position; // mask for the non-representative cluster.
		this->read_cluster_bit_representation ^= mask2; // XOR operation to flip the bit.
	} else if (graycode_bit_changed != -1) {
		// Updating b_index using masks.
		uint32_t mask = 1 << *cluster_bit_changed;
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
	assert(this->bipartition_index < parent->get_num_bipartition());
	return this->bipartition_index;
}

uint32_t BipartitionIterator::get_read_cluster_bit_representation() const {
	return this->read_cluster_bit_representation;
}

bool BipartitionIterator::is_clustered_bit(uint32_t cluster_bit_changed) const {
	return parent->get_cluster_id_to_read_index_map()->count(cluster_bit_changed) > 0;
}

bool BipartitionIterator::has_constrained_bit(uint32_t cluster_bit_changed) const {

	return parent->get_constrained_position_map()->count(cluster_bit_changed) > 0;
}

uint32_t BipartitionIterator::get_constrained_bit(uint32_t cluster_bit_changed) const {
	assert(is_clustered_bit(cluster_bit_changed));
	assert(has_constrained_bit(cluster_bit_changed));
	return parent->get_constrained_position_map()->at(cluster_bit_changed);
}

const std::vector<uint32_t>* BipartitionIterator::get_read_index_from_cluster_id(uint32_t cluster_id) const {
	const std::unordered_map<uint32_t, std::vector<uint32_t>> cluster_id_to_read_index_map = *(parent->get_cluster_id_to_read_index_map());
	return &(cluster_id_to_read_index_map.at(cluster_id));
}

void BipartitionIterator::get_changed_reads(uint32_t cluster_bit_changed, std::unordered_map<uint32_t, bool>& changed_reads) const {
	bool new_bit = (read_cluster_bit_representation >> cluster_bit_changed) & 1;
	uint32_t flipped_cluster_id = parent->get_read_cluster_ids()->at(cluster_bit_changed);
	if (!is_clustered_bit(cluster_bit_changed)) {
		// The bit did not correspond to a cluster.
		changed_reads[flipped_cluster_id] = new_bit; // cluster ID is the read ID.
		return ;
	}
	/**
	 * The bit changed corresponds to a cluster.
	 * Adding the information of the read_ids corresponding to this cluster.
	 */
	const std::vector<uint32_t>& read_indices_from_clusters = *get_read_index_from_cluster_id(flipped_cluster_id);
	for (auto read_index: read_indices_from_clusters) {
		changed_reads[read_index] = new_bit;
	}
	/**
	 * If there is a constrained cluster, add the reads id from that cluster.
	 */
	if (has_constrained_bit(cluster_bit_changed)) {
		uint32_t constrained_flipped_cluster_id = parent->get_read_cluster_ids()->at(parent->get_constrained_position_map()->at(cluster_bit_changed));
		std::vector<uint32_t> constrained_read_index_from_clusters = *get_read_index_from_cluster_id(constrained_flipped_cluster_id);
		for (auto read_index: constrained_read_index_from_clusters) {
			changed_reads[read_index] = !new_bit;
		}
	}
}