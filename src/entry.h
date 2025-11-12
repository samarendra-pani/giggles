// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#ifndef ENTRY_H
#define ENTRY_H

#include <iostream>
#include<vector>

class Entry {
	
	public:
		typedef enum { REF_ALLELE = 0, ALT_ALLELE = 1, BLANK = 2, EQUAL_SCORES = 3 } allele_t;
		
		Entry(uint32_t r, uint32_t a, std::vector<uint32_t> s);
		Entry(uint32_t r, allele_t a, std::vector<uint32_t> s);
		Entry();
		
		void set_read_id(uint32_t r) { read_id = r; }
		void set_allele(uint32_t a) { allele = a; }
		void set_scores(const std::vector<uint32_t>& s) { scores = s; }
		
		uint32_t get_read_id() const { return read_id; }
		uint32_t get_allele() const { return allele; }
		std::vector<uint32_t> get_scores() const { return scores; }
		uint32_t get_phred_score() const { return 30; } // Currently hardcoded to 30 since Whatshap uses fixed quality scores for phasing.

		allele_t get_allele_type() const;
		void convert_scores_to_probability();
		//void convert_scores_to_softmin_probability(uint32_t temperature);

		friend std::ostream& operator<<(std::ostream& out, const Entry& e);

	private:
		bool is_sv; // is this a structural variant? If not, then this will be used for the phasing algorithm.
		uint32_t read_id; // zero-based read identifier
		uint32_t allele; // allele type
		std::vector<uint32_t> scores; // distance scores for all alleles
		std::vector<long double> emission_scores; // emission probabilities for all alleles used for genotyping
};

#endif
