/*
Taken from whatshap (version 2.8)
Original filename: src/pedigreedcolumncostcomputer.h
*/

#ifndef PHASING_COLUMN_COST_COMPUTER_H
#define PHASING_COLUMN_COST_COMPUTER_H

#include <array>
#include <vector>
#include <set>
#include <memory>
#include <map>
#include <utility>
#include <array>

#include "../entry.h"
#include "phasingcolumnindexingiterator.h"
#include "../genotypelikelihoods.h"
#include "../genotypingalgorithm.h"


class PhasingColumnCostComputer {
	
	private:
		const std::vector<const Entry*>& column;
		size_t column_index;
		unsigned int partitioning;
		std::vector<std::array<unsigned int, 2>> cost_partition;
		typedef struct allele_assignment_t {
			/** The i-th bit in assignment gives the allele assigned to pedigree partition i. */
			unsigned int assignment;
			/** Cost of this assignment incurred by genotype changes. */
			unsigned int cost;
			allele_assignment_t() : assignment(0), cost(0) {}
			allele_assignment_t(unsigned int assignment, unsigned int cost) : assignment(assignment), cost(cost) {}
		} allele_assignment_t;
		/** All allowed assignments and their costs. */
		std::vector<allele_assignment_t> allele_assignments;
		const std::vector<GenotypingAlgorithm::variant_information_t>* variant_info_table;

	public:

		PhasingColumnCostComputer(const std::vector<const Entry*>& column, size_t column_index, const std::vector<GenotypingAlgorithm::variant_information_t>* variant_info_table);

		void set_partitioning(unsigned int partitioning);

		void update_partitioning(int bit_to_flip);

		unsigned int get_cost();

		typedef struct phased_variant_t {
			Entry::allele_t allele0;
			Entry::allele_t allele1;
			unsigned int quality;
			phased_variant_t() : allele0(Entry::BLANK), allele1(Entry::BLANK), quality(0) {}
			phased_variant_t(Entry::allele_t allele0, Entry::allele_t allele1) : allele0(allele0), allele1(allele1), quality(0) {}
			phased_variant_t(Entry::allele_t allele0, Entry::allele_t allele1, unsigned int quality) : allele0(allele0), allele1(allele1), quality(quality) {}
		} phased_variant_t;

		/** Returns a phased variants. */
		phased_variant_t get_alleles();

};

#endif
