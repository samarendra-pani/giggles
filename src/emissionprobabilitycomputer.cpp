#include <cassert>

#include "emissionprobabilitycomputer.h"

EmissionProbabilityComputer::EmissionProbabilityComputer(uint32_t n_alleles): i_multiplier(n_alleles, 1.0L), j_multiplier(n_alleles, 1.0L) {
    emission_probability_table = Vector2D<long double>(n_alleles, n_alleles, 1.0L);
}

long double EmissionProbabilityComputer::at(uint32_t i, uint32_t j) const {
    return emission_probability_table.at(i, j);
}

void EmissionProbabilityComputer::update_emission_probability(const int cluster_bit_changed, const BipartitionIterator& iterator, std::vector<const Entry *>& entries) {
    uint32_t n_alleles = emission_probability_table.get_size0();
	long double* table_data = emission_probability_table.data();
	if (cluster_bit_changed >= 0) {
		/**
		 * A cluster has been flipped since the last call to this function.
		 */
		
		/**
		 * i_multiplier represents the 0-th partition and j_multiplier represents the 1-st partition.
		 * If new bit is 1, then that means the partition changed from 0 to 1.
		 * 		So j_multiplier gets the emissions and i_multiplier gets the reciprocal emissions.
		 * If new bit is 0, then that means the partition changed from 1 to 0.
		 * 		So i_multiplier gets the emissions and j_multiplier gets the reciprocal emissions.
		 */
		std::fill(i_multiplier.begin(), i_multiplier.begin() + n_alleles, 1.0L);
    	std::fill(j_multiplier.begin(), j_multiplier.begin() + n_alleles, 1.0L);

		/* Determine numerator and denominator multipliers across ALL changed reads */
		changed_reads.clear();
		iterator.get_changed_reads(cluster_bit_changed, changed_reads);
		for (auto const& [entry_index, newBit]: changed_reads) {
			const Entry* entry = entries[entry_index];
			if (entry->get_allele_type() == Entry::BLANK) continue;
			if (newBit) {
				for (uint32_t allele = 0; allele < n_alleles; allele++) {
					i_multiplier[allele] *= entry->get_reciprocal_emission_score(allele);
					j_multiplier[allele] *= entry->get_emission_score(allele);
				}
			}
			else {
				for (uint32_t allele = 0; allele < n_alleles; allele++) {
					i_multiplier[allele] *= entry->get_emission_score(allele);
					j_multiplier[allele] *= entry->get_reciprocal_emission_score(allele);
				}
			}
		}

		/* Update emissions using the multipliers. */
		for (uint32_t i = 0; i < n_alleles; i++) {
			long double i_mult = i_multiplier[i];
			for (uint32_t j = 0; j < n_alleles; j++) {
				long double total_ratio = i_mult * j_multiplier[j];
				uint32_t flat_index = i*n_alleles + j;
				table_data[flat_index] = table_data[flat_index] * total_ratio;
			}
		}
	}
	else {
		/**
		 * Initialization case.
		 * The  is at the first bipartition.
		 */
		
		/* Separate the reads into bipartition 0 and 1 */
		std::vector<uint32_t> bipar_0_reads;
		std::vector<uint32_t> bipar_1_reads;

		uint32_t cluster_count = 0;
		uint32_t read_cluster_bit_representation = iterator.get_read_cluster_bit_representation();
		const std::vector<uint32_t>* cluster_ids = iterator.get_parent_column()->get_read_cluster_ids();

		for (uint32_t const c_id: *cluster_ids) {
			bool bit = (read_cluster_bit_representation & (1 << cluster_count)) != 0;
			const std::vector<uint32_t>& read_indices = iterator.get_read_index_from_cluster_id(c_id);
			for (uint32_t const r_idx: read_indices) {
				if (entries[r_idx]->get_allele_type() != Entry::BLANK) {
					if (bit) bipar_1_reads.push_back(r_idx);
					else     bipar_0_reads.push_back(r_idx);
				}
			}
			cluster_count++;
		}

		/* Compute independent products for bipartition 0 and 1 */
		std::vector<long double> prod_0(n_alleles, 1.0L);
		std::vector<long double> prod_1(n_alleles, 1.0L);
		for (uint32_t i = 0; i < n_alleles; i++) {
			for (uint32_t r_idx : bipar_0_reads) {
				prod_0[i] *= entries[r_idx]->get_emission_score(i);
			}
			for (uint32_t r_idx : bipar_1_reads) {
				prod_1[i] *= entries[r_idx]->get_emission_score(i);
			}
		}

		/* Calculate emission probabilities */
		for (uint32_t i = 0; i < n_alleles; i++) {
			for (uint32_t j = 0; j < n_alleles; j++) {
				// emission_probability_table.set(i, j, prod_0[i] * prod_1[j]);
				uint32_t flat_index = i*n_alleles + j;
				table_data[flat_index] = prod_0[i]*prod_1[j];
			}
		}
	}
}
    