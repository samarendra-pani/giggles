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

PhasingColumnCostComputer::PhasingColumnCostComputer(const std::vector <const Entry *>&column, size_t column_index, const std::vector<GenotypingAlgorithm::variant_information_t>* variant_info_table):
	column(column),
	column_index(column_index),
	variant_info_table(variant_info_table),
	partitioning(0)
{
	// Enumerate all possible assignments of alleles to haplotypes and 
	// store those that are compatible with genotypes.
	for (unsigned int i = 0; i < (1<<2); ++i) {
		bool genotypes_compatible = true;
		unsigned int cost = 0;
		unsigned int allele0 = (i >> 0) & 1;
		unsigned int allele1 = (i >> 1) & 1;
		Genotype genotype(vector<unsigned int>{allele0,allele1});
		const GenotypeLikelihoods* gls = &variant_info_table->at(column_index).genotype_likelihoods;
		assert(gls != nullptr);
		cost += gls->getPhredScore(genotype);
		if (genotypes_compatible) {
			allele_assignments.push_back(allele_assignment_t(i,cost));
		}
	}
}


void PhasingColumnCostComputer::set_partitioning(unsigned int partitioning) {
	cost_partition.assign(2, {0,0});	// two partitions, each with cost for ref and alt

	partitioning = partitioning;
	for (vector < const Entry * >::const_iterator it = column.begin(); it != column.end(); ++it) {
		auto & entry = **it;
		bool entry_in_partition1 = (partitioning & ((unsigned int) 1)) == 0;
		switch (entry.get_allele_type()) {
			case Entry::REF_ALLELE:
				(entry_in_partition1 ? cost_partition[0] :cost_partition[1])[1] += entry.get_phred_score();
				break;
			case Entry::ALT_ALLELE:
				(entry_in_partition1 ? cost_partition[0] :cost_partition[1])[0] += entry.get_phred_score();
				break;
			case Entry::BLANK:
				break;
			default:
				assert(false);
		}
		partitioning = partitioning >> 1;
	}
}


void PhasingColumnCostComputer::update_partitioning(int bit_to_flip) {
	const Entry & entry = *column[bit_to_flip];
	partitioning = partitioning ^ (((unsigned int) 1) << bit_to_flip);
	bool entry_in_partition1 = (partitioning & (((unsigned int) 1) << bit_to_flip)) == 0;
	unsigned int ind_id = 0; // only one individual in the pedigree
	switch (entry.get_allele_type()) {
		case Entry::REF_ALLELE:
			(entry_in_partition1 ? cost_partition[1] : cost_partition[0])[1] -= entry.get_phred_score();
			(entry_in_partition1 ? cost_partition[0] :  cost_partition[1])[1] += entry.get_phred_score();
			break;
		case Entry::ALT_ALLELE:
			(entry_in_partition1 ? cost_partition[1] : cost_partition[0])[0] -= entry.get_phred_score();
			(entry_in_partition1 ? cost_partition[0] :  cost_partition[1])[0] += entry.get_phred_score();
			break;
		case Entry::BLANK:
			break;
		default:
			assert(false);
	}
}


unsigned int PhasingColumnCostComputer::get_cost() {
	unsigned int best_cost = numeric_limits < unsigned int >::max();
	for (const allele_assignment_t& a : allele_assignments) {
		unsigned int cost = a.cost;
		// there are only two partitions
		for (size_t p = 0; p < 2; ++p) {
			unsigned int allele = (a.assignment >> p) & 1;
			cost += cost_partition[p][allele];
		}
		if (cost < best_cost) {
			best_cost = cost;
		}
	}
	return best_cost;
}


PhasingColumnCostComputer::phased_variant_t PhasingColumnCostComputer::get_alleles() {
	unsigned int best_cost = numeric_limits < unsigned int >::max();
	unsigned int second_best_cost = numeric_limits < unsigned int >::max();
	phased_variant_t haps;
	// best_cost_for_allele[haplotype][to_allele] is the best cost for flipping
	vector<array<unsigned int,2>> best_cost_for_allele(2, {numeric_limits<unsigned int>::max(),numeric_limits<unsigned int>::max()});
	for (const allele_assignment_t& a : allele_assignments) {
		unsigned int cost = a.cost;
		for (size_t p = 0; p < 2; ++p) {
			unsigned int allele = (a.assignment >> p) & 1;
			cost += cost_partition[p][allele];
		}
		bool new_best = false;
		if (cost <= best_cost) {
			best_cost = cost;
			new_best = true;
		}
		unsigned int partition0 = 0;
		unsigned int partition1 = 1;
		unsigned int allele0 = (a.assignment >> 0) & 1;
		unsigned int allele1 = (a.assignment >> partition1) & 1;
		if (new_best) {
			haps = phased_variant_t(
				(allele0 == 0)?Entry::REF_ALLELE:Entry::ALT_ALLELE,
				(allele1 == 0)?Entry::REF_ALLELE:Entry::ALT_ALLELE
			);
		}
		if (cost < best_cost_for_allele.at(0)[allele0]) {
			best_cost_for_allele.at(0)[allele0] = cost;
		}
		if (cost < best_cost_for_allele.at(1)[allele1]) {
			best_cost_for_allele.at(1)[allele1] = cost;
		}
		
	}

	if (best_cost == numeric_limits < unsigned int >::max()) {
		throw std::runtime_error("Error: Mendelian conflict");
	}

	// Test whether some of the allele assignments are ambiguous
	for (size_t haplotype = 0; haplotype < 2; ++haplotype) {
		int quality = abs(((int)(best_cost_for_allele.at(haplotype)[0])) - ((int)(best_cost_for_allele.at(haplotype)[1])));
		haps.quality = (unsigned int)quality;
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