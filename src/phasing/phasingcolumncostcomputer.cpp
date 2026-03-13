/*
Taken from whatshap (version 2.8)
Original filename: src/pedigreedcolumncostcomputer.cpp
*/

#include <cassert>
#include <limits>
#include <utility>
#include <algorithm>
#include <array>
#include <map>
#include <unordered_set>
#include "../vector2d.h"
#include "../genotype.h"
#include "phasingcolumncostcomputer.h"

using namespace std;

PhasingColumnCostComputer::PhasingColumnCostComputer(const std::vector <const Entry *>&column, const variant_information_t& variant_info):
	column(column),
	partitioning(0) {

	/** Ensures that at most 2 alleles are active. */
	if (!variant_info.phasable) {
		this->phasable = false;
		return;
	}

	/** Pileup analysis of the entry alleles to get some idea about which genotype it might be. */
	std::vector<uint32_t> selected_genotype_indices = variant_info.genotype_likelihoods.select_genotypes();
	std::vector<uint32_t> compatible_genotype_indices;
	assert(selected_genotype_indices.size() <= 3); // With <= 2 active alleles, there can only be a max of 3 selected genotypes possible.
	if (selected_genotype_indices.size() == 3) {
		/** All genotype combinations seems to be possible. Doing a pileup analysis of entries */
		switch (analyse_entry_alleles()) {
			case 0:
				compatible_genotype_indices.push_back(selected_genotype_indices[0]); // Adding index for HOM-ALLELE1
				break;
			case 1:
				compatible_genotype_indices.push_back(selected_genotype_indices[0]); // Adding index for HOM-ALLELE1
				compatible_genotype_indices.push_back(selected_genotype_indices[1]); // Adding index for HET
				break;
			case 2:
				compatible_genotype_indices.push_back(selected_genotype_indices[1]); // Adding index for HET
				break;
			case 3:
				compatible_genotype_indices.push_back(selected_genotype_indices[1]); // Adding index for HET
				compatible_genotype_indices.push_back(selected_genotype_indices[2]); // Adding index for HOM-ALLELE2
				break;
			case 4:
				compatible_genotype_indices.push_back(selected_genotype_indices[2]); // Adding index for HOM-ALLELE2
				break;
			case -1:
				compatible_genotype_indices = selected_genotype_indices; // Adding all genotype indices.
				break;
			default:
				// there should not be any other output
				assert(false);
				break;
		}
	} else {
		/** Adding the genotype indices selected through genotype likelihoods. */
		compatible_genotype_indices = selected_genotype_indices;
	}
	/**
	 * Enumerate all possible assignments of alleles to haplotypes and 
	 * store those that are compatible with genotypes.
	 * 
	 * Let's say we have Haplotypes 0 and 1 and each Haplotype can be either Allele 0 or 1 (where these two are the active alleles.)
	 * 
	 * i = 0 -> Hap0 has Allele0 (allele0 = 0) and Hap1 has Allele0 (allele1 = 0)
	 * i = 1 -> Hap0 has Allele1 (allele0 = 1) and Hap1 has Allele0 (allele1 = 0)
	 * i = 2 -> Hap0 has Allele0 (allele0 = 0) and Hap1 has Allele1 (allele1 = 1)
	 * i = 3 -> Hap0 has Allele1 (allele0 = 1) and Hap1 has Allele1 (allele1 = 1)
	 */
	std::vector<uint32_t> active_alleles = variant_info.get_active_positions();
	uint32_t allele0;
	uint32_t allele1;
	Genotype genotype;
	uint32_t index;
	bool genotypes_compatible;
	if (active_alleles.size() == 1) {
		assert(compatible_genotype_indices.size() == 1); // There should be only compatible genotype -> homozygous of the active allele.
		allele_assignments.push_back(0);
	} else {
		assert(active_alleles.size() == 2);
		for (uint32_t i = 0; i < 4; ++i) {
			allele0 = active_alleles[(i >> 0) & 1];
			allele1 = active_alleles[(i >> 1) & 1];
			genotype = Genotype(vector<uint32_t>{allele0,allele1});
			index = genotype.get_index();
			/** Checking if index exists in the compatible genotype indices. */
			if (std::find(compatible_genotype_indices.begin(), compatible_genotype_indices.end(), index) != compatible_genotype_indices.end()) {
				allele_assignments.push_back(i);
			}
		}
	}
}


void PhasingColumnCostComputer::set_partitioning(uint32_t partitioning) {
	if (!phasable) {
		return;
	}

	cost_partition.assign(2, {0,0});	// two partitions, each with cost for ref and alt

	partitioning = partitioning;
	for (vector < const Entry * >::const_iterator it = column.begin(); it != column.end(); ++it) {
		auto & entry = **it;
		bool entry_in_partition1 = (partitioning & ((uint32_t) 1)) == 0;
		switch (entry.get_allele_type()) {
			case Entry::ALLELE1:
				(entry_in_partition1 ? cost_partition[0] :cost_partition[1])[1] += entry.get_phred_score();
				break;
			case Entry::ALLELE2:
				(entry_in_partition1 ? cost_partition[0] :cost_partition[1])[0] += entry.get_phred_score();
				break;
			case Entry::BLANK:
				break;
			case Entry::EQUAL_SCORES:
				break;
			default:
				assert(false);
		}
		partitioning = partitioning >> 1;
	}
}


void PhasingColumnCostComputer::update_partitioning(int bit_to_flip) {
	if (!phasable) {
		return;
	}
	const Entry & entry = *column[bit_to_flip];
	partitioning = partitioning ^ (((uint32_t) 1) << bit_to_flip);
	bool entry_in_partition1 = (partitioning & (((uint32_t) 1) << bit_to_flip)) == 0;
	uint32_t ind_id = 0; // only one individual in the pedigree
	switch (entry.get_allele_type()) {
		case Entry::ALLELE1:
			(entry_in_partition1 ? cost_partition[1] : cost_partition[0])[1] -= entry.get_phred_score();
			(entry_in_partition1 ? cost_partition[0] :  cost_partition[1])[1] += entry.get_phred_score();
			break;
		case Entry::ALLELE2:
			(entry_in_partition1 ? cost_partition[1] : cost_partition[0])[0] -= entry.get_phred_score();
			(entry_in_partition1 ? cost_partition[0] :  cost_partition[1])[0] += entry.get_phred_score();
			break;
		case Entry::BLANK:
			break;
		case Entry::EQUAL_SCORES:
			break;
		default:
			assert(false);
	}
}


uint32_t PhasingColumnCostComputer::get_cost() {
	if (!phasable) {
		return 0;
	}
	uint32_t best_cost = numeric_limits < uint32_t >::max();
	for (const uint32_t a : allele_assignments) {
		uint32_t cost = 0;
		// there are only two partitions
		for (size_t p = 0; p < 2; ++p) {
			uint32_t allele = (a >> p) & 1;
			cost += cost_partition[p][allele];
		}
		if (cost < best_cost) {
			best_cost = cost;
		}
	}
	return best_cost;
}


PhasingColumnCostComputer::phased_variant_t PhasingColumnCostComputer::get_alleles() {
	
	phased_variant_t haps;
	if (!phasable) {
		haps.allele0 = Entry::EQUAL_SCORES;
		haps.allele1 = Entry::EQUAL_SCORES;
		haps.quality = 0;
		return haps;
	}
	uint32_t best_cost = numeric_limits < uint32_t >::max();
	uint32_t second_best_cost = numeric_limits < uint32_t >::max();
	
	// best_cost_for_allele[haplotype][to_allele] is the best cost for flipping
	vector<array<uint32_t,2>> best_cost_for_allele(2, {numeric_limits<uint32_t>::max(),numeric_limits<uint32_t>::max()});
	for (const uint32_t a : allele_assignments) {
		uint32_t cost = 0;
		for (size_t p = 0; p < 2; ++p) {
			uint32_t allele = (a >> p) & 1;
			cost += cost_partition[p][allele];
		}
		bool new_best = false;
		if (cost <= best_cost) {
			best_cost = cost;
			new_best = true;
		}
		uint32_t partition0 = 0;
		uint32_t partition1 = 1;
		uint32_t allele0 = (a >> 0) & 1;
		uint32_t allele1 = (a >> partition1) & 1;
		if (new_best) {
			haps = phased_variant_t(
				(allele0 == 0)?Entry::ALLELE1:Entry::ALLELE2,
				(allele1 == 0)?Entry::ALLELE1:Entry::ALLELE2
			);
		}
		if (cost < best_cost_for_allele.at(0)[allele0]) {
			best_cost_for_allele.at(0)[allele0] = cost;
		}
		if (cost < best_cost_for_allele.at(1)[allele1]) {
			best_cost_for_allele.at(1)[allele1] = cost;
		}
		
	}

	if (best_cost == numeric_limits < uint32_t >::max()) {
		throw std::runtime_error("Error: Mendelian conflict");
	}

	// Test whether some of the allele assignments are ambiguous
	for (size_t haplotype = 0; haplotype < 2; ++haplotype) {
		uint32_t quality = abs(((int)(best_cost_for_allele.at(haplotype)[0])) - ((int)(best_cost_for_allele.at(haplotype)[1])));
		haps.quality = (uint32_t)quality;
		if (quality == 0) {
			if (haplotype == 0) {
				haps.allele0 = Entry::EQUAL_SCORES;
			} else {
				haps.allele1 = Entry::EQUAL_SCORES;
			}
		}
	}

	return haps;
}

int PhasingColumnCostComputer::analyse_entry_alleles() {
	uint32_t total_count = column.size();
	// if not enough entries
	if (total_count < 5) {
		return -1;
	}
	uint32_t allele1_count = 0;
	uint32_t allele2_count = 0;
	uint32_t unknown_count = 0;
	for (const Entry* e: column) {
		switch (e->get_allele_type()) {
			case Entry::ALLELE1:
				allele1_count++;
				break;
			case Entry::ALLELE2:
				allele2_count++;
				break;
			default:
				unknown_count++;
				break;
		}
	}
	// ratio of informative enties is low
	if ((allele1_count + allele2_count)/total_count < 0.7) {
		return -1;
	}
	double allele1_frac = (allele1_count)/((allele1_count + allele2_count));
	if (allele1_frac > 0.9) { return 0; } // HOM-ALLELE1
	if (allele1_frac > 0.6) { return 1; } // HET/HOM-ALLELE1
	if (allele1_frac > 0.4) { return 2; } // HET
	if (allele1_frac > 0.1) { return 3; }  // HET/HOM-ALLELE2
	return 4; // HOM-ALLELE2
}