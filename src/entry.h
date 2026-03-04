// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

/*
Since WhatsHap phasing is restricted to two alleles, we have to hard code ALLELE1 and ALLELE2.
*/

#ifndef ENTRY_H
#define ENTRY_H

#include <iostream>
#include <vector>

#include "genotype.h"

class Entry {
	
	public:
		typedef enum: uint8_t { ALLELE1 = 0, ALLELE2 = 1, BLANK = 2, EQUAL_SCORES = 3 } allele_t;
		
		Entry(uint32_t r, const std::vector<uint32_t>& s);
		Entry();
		
		void set_read_id(uint32_t r);
		void set_scores(const std::vector<uint32_t>& s);
		/* 
		* set allele type based on the active alleles at that position.
		*   if the allele at the lower index value has lower distance score, then assigned ALLELE1.
		*   if the allele at the higher index value has lower distance score, then assigned ALLELE2.
		*   if both have same distance score, then assigned EQUAL_SCORES.
		*/
		void set_allele_type(const std::vector<bool>& active_alleles);
		/*
		* set allele type based on pre-computed allele type.
		* used for the super-reads created to represent the haplotypes.
		*/
		void set_allele_type(allele_t a);
		
		uint32_t get_read_id() const;
		uint32_t get_phred_score() const;
		/**
		 * Return allele type: ALLELE1, ALLELE2, BLANK, or EQUAL_SCORES
		 */
		allele_t get_allele_type() const;
		/**
		 * Returns false if allele_type is BLANK
		 */
		bool has_allele_type() const;
		
		void convert_scores_to_probability(const std::vector<uint32_t>& scores);
		//void convert_scores_to_softmin_probability(uint32_t temperature);

		std::vector<long double> get_emission_scores() const;
		void set_emission_scores(const std::vector<long double>& scores);

		friend std::ostream& operator<<(std::ostream& out, const Entry& e);

	private:
		uint32_t read_id;	// zero-based read identifier
		allele_t allele;	// allele type
		std::vector<long double> emission_scores; // emission probabilities for all alleles used for genotyping
};

#endif
