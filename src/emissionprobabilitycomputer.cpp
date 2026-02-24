#include <cassert>

#include "emissionprobabilitycomputer.h"

EmissionProbabilityComputer::EmissionProbabilityComputer(uint32_t n_alleles) {
    emission_probability_table = Vector2D<long double>(n_alleles, n_alleles, 1.0L);
}

long double EmissionProbabilityComputer::at(uint32_t i, uint32_t j) const {
    return emission_probability_table.at(i, j);
}

void EmissionProbabilityComputer::update_emission_probability(const int cluster_bit_changed, const BipartitionIterator& iterator, std::vector<const Entry *>& entries) {
    uint32_t n_alleles = emission_probability_table.get_size0();
	if (cluster_bit_changed >= 0) {
		/**
		 * A cluster has been flipped since the last call to this function.
		 */
		changed_reads.clear();
		iterator.get_changed_reads(cluster_bit_changed, changed_reads);
		long double ratio;
		const Entry* entry;
		for (auto const& [entry_index, newBit]: changed_reads) {
			// Need to check if the entry being changed actually corresponds to the correct read id.
			entry = entries[entry_index];
			assert(entry->get_read_id() == iterator.get_parent_column()->get_read_ids()->at(entry_index));
			if (entry->get_allele_type() == Entry::BLANK) {
				/**
				 * The flipped entry has no effect on emission. It is blank.
				 */
				continue;
			}
			const std::vector<long double>& scores = entries.at(entry_index)->get_emission_scores();
			for (uint32_t i = 0; i < n_alleles; i++) {
				for (uint32_t j = 0; j < n_alleles; j++) {
					if (newBit) {
						/**
						 * The new bit is 1. So the entry was at bipartion 0 (corresponding to allele i) and is now at bipartition 1 (corresponding to allele j).
						 */
						ratio = scores[j]/scores[i];
					} 
					else {
						/**
						 * The new bit is 0. So the entry was at bipartion 1 (corresponding to allele j) and is now at bipartition 0 (corresponding to allele i).
						 */
						ratio = scores[i]/scores[j];
					}
					emission_probability_table.set(i, j, emission_probability_table.at(i, j) * ratio);
				}
			}
		}	
	}
	else {
		/**
		 * Initialization case.
		 * The  is at the first bipartition.
		 */
		uint32_t read_cluster_bit_representation = iterator.get_read_cluster_bit_representation(); // this is value showing the bipartition of every read cluster.
		/**
		 * Iterating through all allele pairs
		 */
		uint32_t cluster_count;
		long double value;
		const std::vector<uint32_t>* cluster_ids;
		bool bit;
		for (uint32_t i = 0; i < n_alleles; i++) {
			for (uint32_t j = 0; j < n_alleles; j++) {
				cluster_count = 0;
				value = 1.0L;
				/**
				 * Finding the cluster ids in the current column.
				 * The bits of read_cluster_bit_representation correspond to these cluster ids.
				 */
				cluster_ids = iterator.get_parent_column()->get_read_cluster_ids();
				for (uint32_t const c_id: *cluster_ids) {
					/**
					 * For each cluster, get the read indices corresponding to that cluster.
					 * These read indices correspond to the indices of their respective entries.
					 */
					const std::vector<uint32_t>& read_indices = iterator.get_read_index_from_cluster_id(c_id);
					// Determine whether the current cluster is in bipartition 0 or 1
					bit = (read_cluster_bit_representation & (1 << cluster_count)) != 0;
					for (uint32_t const r_idx: read_indices) {
						if (entries[r_idx]->get_allele_type() == Entry::BLANK) {
							continue;
						}
						if (bit) {
							// If read at r_idx is in biparition 1 then value gets multiplied with emission from allele j
							value = value * (entries[r_idx]->get_emission_scores()[j]);
						}
						else {
							// If read at r_idx is in biparition 0 then value gets multiplied with emission from allele i
							value = value * (entries[r_idx]->get_emission_scores()[i]);
						}
					}
					// Move to the next cluster
					cluster_count++;
				}
				// Setting the computed emission probability value
				emission_probability_table.set(i, j, value);
			}
		}
	}
}
    