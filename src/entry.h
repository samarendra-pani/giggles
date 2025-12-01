// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

/*
Since WhatsHap phasing is restricted to two alleles, we have to hard code ALLELE1 and ALLELE2.
*/

#ifndef ENTRY_H
#define ENTRY_H

#include <iostream>
#include <vector>

class Entry {
	
	public:
		typedef enum { ALLELE1 = 0, ALLELE2 = 1, BLANK = 2, EQUAL_SCORES = 3 } allele_t;
		
		Entry(uint32_t r, std::vector<uint32_t> s);
		Entry();
		
		void set_read_id(uint32_t r);
		void set_scores(const std::vector<uint32_t>& s);
		/* 
		* set allele type based on the active alleles at that position.
		*   if the allele at the lower index value has lower distance score, then assigned ALLELE1.
		*   if the allele at the higher index value has lower distance score, then assigned ALLELE2.
		*   if both have same distance score, then assigned EQUAL_SCORES.
		*/
		void set_allele_type(std::vector<bool> active_alleles);
		/*
		* set allele type based on pre-computed allele type.
		* used for the super-reads created to represent the haplotypes.
		* need the information of which active alleles are used for 
		*/
		void set_allele_type(allele_t a, uint32_t idx1, uint32_t idx2);
		
		uint32_t get_read_id() const;
		std::vector<uint32_t> get_scores() const;
		uint32_t get_phred_score() const;
		allele_t get_allele_type() const;

		bool has_allele_type() const;
		
		void convert_scores_to_probability();
		//void convert_scores_to_softmin_probability(uint32_t temperature);

		friend std::ostream& operator<<(std::ostream& out, const Entry& e);

	private:
		uint32_t read_id; // zero-based read identifier
		allele_t allele; // allele type
		uint32_t allele1_idx; // index of allele 1
		uint32_t allele2_idx; // index of allele 2
		std::vector<uint32_t> scores; // distance scores for all alleles
		std::vector<long double> emission_scores; // emission probabilities for all alleles used for genotyping
};

#endif
