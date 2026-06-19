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
		
		Entry(uint32_t r, const std::vector<float>& s);
		Entry(uint32_t r);
		Entry();
		
		void set_read_id(uint32_t r);
		void set_scores(const std::vector<float>& s);
		void set_scores(const std::vector<uint8_t>& s);
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

		long double get_emission_score(uint32_t i) const;
		long double get_reciprocal_emission_score(uint32_t j) const;

		static void initialize_probability_cache(float temperature);

		friend std::ostream& operator<<(std::ostream& out, const Entry& e);

	private:
   		static std::vector<long double> probability_cache;	// cache of probabilities defined by the distance
		static std::vector<long double> reciprocal_probability_cache;	// cache of reciprocal probabilities defined by the distance
		static uint32_t k;	// number of discrete values
		uint32_t read_id;	// zero-based read identifier
		allele_t allele;	// allele type
		/**
		 * scores contains the values from realign() from variants.py
		 * now the values have been renormalized from [0, 1] bound to a [0, 100] bound
		 * and also discretised.
		 */
		std::vector<uint8_t> scores;
};

#endif
