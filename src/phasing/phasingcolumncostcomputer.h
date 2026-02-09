/*
Taken from whatshap (version 2.8)
Original filename: src/pedigreedcolumncostcomputer.h
*/

#ifndef PHASING_COLUMN_COST_COMPUTER_H
#define PHASING_COLUMN_COST_COMPUTER_H

#include <vector>
#include <set>
#include <memory>
#include <map>
#include <utility>
#include <array>

#include "../entry.h"
#include "phasingcolumnindexingiterator.h"
#include "../genotypelikelihoods.h"
#include "../variantinfo.h"


class PhasingColumnCostComputer {
	
	private:
		const std::vector<const Entry*>& column;
		size_t column_index;
		uint32_t partitioning;
		std::vector<std::array<uint32_t, 2>> cost_partition;
		typedef struct allele_assignment_t {
			/** The i-th bit in assignment gives the allele assigned to pedigree partition i. */
			uint32_t assignment;
			/** Cost of this assignment incurred by genotype changes. */
			uint32_t cost;
			allele_assignment_t() : assignment(0), cost(0) {}
			allele_assignment_t(uint32_t assignment, uint32_t cost) : assignment(assignment), cost(cost) {}
		} allele_assignment_t;
		/** All allowed assignments and their costs. */
		std::vector<allele_assignment_t> allele_assignments;
		const std::vector<variant_information_t>* variant_info_table;

	public:

		PhasingColumnCostComputer(const std::vector<const Entry*>& column, size_t column_index, const std::vector<variant_information_t>* variant_info_table);

		void set_partitioning(uint32_t partitioning);

		void update_partitioning(int bit_to_flip);

		uint32_t get_cost();

		typedef struct phased_variant_t {
			Entry::allele_t allele0;
			Entry::allele_t allele1;
			uint32_t quality;
			phased_variant_t() : allele0(Entry::BLANK), allele1(Entry::BLANK), quality(0) {}
			phased_variant_t(Entry::allele_t allele0, Entry::allele_t allele1) : allele0(allele0), allele1(allele1), quality(0) {}
			phased_variant_t(Entry::allele_t allele0, Entry::allele_t allele1, uint32_t quality) : allele0(allele0), allele1(allele1), quality(quality) {}
		} phased_variant_t;

		/** Returns a phased variants. */
		phased_variant_t get_alleles();

};

#endif
