// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#include <cassert>
#include <bits/stdc++.h>
#include <cmath>
using namespace std;

#include "entry.h"

Entry::Entry(unsigned int r, unsigned int a, std::vector<unsigned int> s) : 
	read_id(r), allele(a), scores(s), is_sv(false) {}

Entry::Entry(unsigned int r, allele_t a, std::vector<unsigned int> s) : 
	read_id(r), allele((unsigned int)a), scores(s), is_sv(false) {}

Entry::Entry() : read_id(0), allele(2), scores({}), is_sv(false) {}

/*
Previous conversion from distance scores to emission probabilities
{
	emission_score.resize(e.size());
	int i = 0;
	double normalization = 0.0L;
	for (auto it = e.begin(); it != e.end(); it++, i++) {
		emission_score[i] = pow(10, -reg_const);
		long double score = 0.0L;
		for (auto it2 = e.begin(); it2 != e.end(); it2++) {
			score += pow(base_const, *it2 - *it);
		}
		emission_score[i] += 1/score;
		normalization += emission_score[i];
	}
	transform((emission_score).begin(), (emission_score).end(), (emission_score).begin(), std::bind2nd(std::divides<long double>(), normalization));
}
*/	

Entry::allele_t Entry::get_allele_type() const {
	if (scores.size() > 1) {
		// this record is multi-allelic, so we cannot assign a single allele type
		// this should never happen in the phasing algorithm
		throw std::runtime_error("Error: cannot determine allele type for multi-allelic variant.");
	}
	switch (allele) {
		case 0: return REF_ALLELE;
		case 1: return ALT_ALLELE;
		case 2: return BLANK;
		case 3: return EQUAL_SCORES;
		default: throw std::runtime_error("Error: invalid allele type.");
	}
}

std::ostream& operator<<(std::ostream& out, const Entry& e) {
	out << "Entry(" << e.read_id ;
	out << ","<< e.allele << ",(";
	for (auto i : e.scores) {
		out << i << ",";
	}
	out << ")" << std::endl;
	return out;
}

