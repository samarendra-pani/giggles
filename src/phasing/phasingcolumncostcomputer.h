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
		uint32_t partitioning;
		bool phasable;
		bool homozygous;
		Entry::allele_t homozygous_allele;
		std::vector<std::array<uint32_t, 2>> cost_partition;
		
		/** All allowed assignments. */
		std::vector<uint32_t> allele_assignments;
		/** 
		 * A pile-up analysis of Entries to check if there is enough evidence to support HOM or HET.
		 * Returns 
		 * 	0 for HOM-ALLELE1
		 * 	1 for HET/HOM-ALLELE1
		 * 	2 for HET
		 * 	3 for HET/HOM-ALLELE2
		 * 	4 for HOM-ALLELE2
		 * 	-1 for insufficient evidence.
		 * 	-2 for absolutely no evidence (all are BLANK or EQUAL_SCORES)
		 * 
		 */
		int analyse_entry_alleles();

	public:

		PhasingColumnCostComputer(const std::vector<const Entry*>& column, const variant_information_t& variant_info);

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
