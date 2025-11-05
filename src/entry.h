// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#ifndef ENTRY_H
#define ENTRY_H

#include <iostream>
#include<vector>

class Entry {
	
	public:
		typedef enum { REF_ALLELE = 0, ALT_ALLELE = 1, BLANK = 2, EQUAL_SCORES = 3 } allele_t;
		
		Entry(unsigned int r, unsigned int a, std::vector<unsigned int> s);
		Entry(unsigned int r, allele_t a, std::vector<unsigned int> s);
		Entry();
		
		void set_read_id(unsigned int r) { read_id = r; }
		void set_allele(unsigned int a) { allele = a; }
		void set_scores(const std::vector<unsigned int>& s) { scores = s; }
		
		unsigned int get_read_id() const { return read_id; }
		unsigned int get_allele() const { return allele; }
		std::vector<unsigned int> get_scores() const { return scores; }
		unsigned int get_phred_score() const { return 30; } // Currently hardcoded to 30 since Whatshap uses fixed quality scores for phasing.

		allele_t get_allele_type() const;
		void convert_scores_to_probability();
		//void convert_scores_to_softmin_probability(unsigned int temperature);

		friend std::ostream& operator<<(std::ostream& out, const Entry& e);

	private:
		bool is_sv; // is this a structural variant? If not, then this will be used for the phasing algorithm.
		unsigned int read_id; // zero-based read identifier
		unsigned int allele; // allele type
		std::vector<unsigned int> scores; // distance scores for all alleles
		std::vector<long double> emission_scores; // emission probabilities for all alleles used for genotyping
};

#endif
