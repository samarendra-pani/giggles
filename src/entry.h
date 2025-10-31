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
		void set_is_sv(bool is_sv) { this->is_sv = is_sv; }
		void set_phred_scores(const std::vector<unsigned int>& s) { phred_scores = s; }
		
		unsigned int get_read_id() const { return read_id; }
		unsigned int get_allele() const { return allele; }
		std::vector<unsigned int> get_scores() const { return scores; }
		bool get_is_sv() const { return is_sv; }
		unsigned int get_phred_score() const { return 30; } // Currently hardcoded to 30 since Whatshap uses fixed quality scores for phasing.

		allele_t get_allele_type() const;

		friend std::ostream& operator<<(std::ostream& out, const Entry& e);

	private:
		bool is_sv; // is this a structural variant? If not, then this will be used for the phasing algorithm.
		unsigned int read_id; // zero-based read identifier
		unsigned int allele; // allele type
		std::vector<unsigned int> scores; // distance scores for all alleles
		std::vector<unsigned int> phred_scores; // phred-scaled scores for all alleles
};

#endif
