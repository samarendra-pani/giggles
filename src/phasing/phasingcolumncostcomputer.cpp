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

PhasingColumnCostComputer::PhasingColumnCostComputer(const std::vector <const Entry *>&column, size_t column_index, const std::vector<variant_information_t>* variant_info_table):
	column(column),
	column_index(column_index),
	variant_info_table(variant_info_table),
	partitioning(0)
{
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
	std::vector<uint32_t> active_alleles = variant_info_table->at(column_index).get_active_positions();
	assert (active_alleles.size() ==  2); // Only two should be active.
	for (uint32_t i = 0; i < 4; ++i) {
		bool genotypes_compatible = true;
		uint32_t cost = 0;
		uint32_t allele0 = active_alleles[(i >> 0) & 1];
		uint32_t allele1 = active_alleles[(i >> 1) & 1];
		Genotype genotype(vector<uint32_t>{allele0,allele1});
		const GenotypeLikelihoods& gls = variant_info_table->at(column_index).genotype_likelihoods;
		assert(gls.size() != 0);
		cost += gls.getPhredScore(genotype);
		allele_assignments.push_back(allele_assignment_t(i,cost));
	}
}


void PhasingColumnCostComputer::set_partitioning(uint32_t partitioning) {
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
	uint32_t best_cost = numeric_limits < uint32_t >::max();
	for (const allele_assignment_t& a : allele_assignments) {
		uint32_t cost = a.cost;
		// there are only two partitions
		for (size_t p = 0; p < 2; ++p) {
			uint32_t allele = (a.assignment >> p) & 1;
			cost += cost_partition[p][allele];
		}
		if (cost < best_cost) {
			best_cost = cost;
		}
	}
	return best_cost;
}


PhasingColumnCostComputer::phased_variant_t PhasingColumnCostComputer::get_alleles() {
	uint32_t best_cost = numeric_limits < uint32_t >::max();
	uint32_t second_best_cost = numeric_limits < uint32_t >::max();
	phased_variant_t haps;
	// best_cost_for_allele[haplotype][to_allele] is the best cost for flipping
	vector<array<uint32_t,2>> best_cost_for_allele(2, {numeric_limits<uint32_t>::max(),numeric_limits<uint32_t>::max()});
	for (const allele_assignment_t& a : allele_assignments) {
		uint32_t cost = a.cost;
		for (size_t p = 0; p < 2; ++p) {
			uint32_t allele = (a.assignment >> p) & 1;
			cost += cost_partition[p][allele];
		}
		bool new_best = false;
		if (cost <= best_cost) {
			best_cost = cost;
			new_best = true;
		}
		uint32_t partition0 = 0;
		uint32_t partition1 = 1;
		uint32_t allele0 = (a.assignment >> 0) & 1;
		uint32_t allele1 = (a.assignment >> partition1) & 1;
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