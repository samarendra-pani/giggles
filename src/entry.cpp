// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#include <cassert>
#include <bits/stdc++.h>
#include <cmath>
using namespace std;

#include "entry.h"

Entry::Entry(unsigned int r, unsigned int a, std::vector<unsigned int> s) : 
	read_id(r), allele(a), scores(s), is_sv(false) {
		convert_scores_to_probability();
	}

Entry::Entry(unsigned int r, allele_t a, std::vector<unsigned int> s) : 
	read_id(r), allele((unsigned int)a), scores(s), is_sv(false) {
		convert_scores_to_probability();
	}

Entry::Entry() : read_id(0), allele(2), scores({}), is_sv(false) {}

/*
* conversion of distance scores to emission probabilities using error probability 0.0001
* emission probability = 10^(-max(score*log10(0.0001), 1e-10))
* modelling the probability of observing a read given the true allele and error rate
*/
void Entry::convert_scores_to_probability() {
	assert(scores.size() > 0);
	long double sum_scores = 0.0L;
	emission_scores.resize(scores.size());
	for (size_t i = 0; i < scores.size(); i++) {
		long double logprob = (long double)std::max(scores[i]*log10(0.0001), 1e-10);
		emission_scores[i] = pow(10.0L, -logprob);
		sum_scores += emission_scores[i];
	}
	// normalizing the emission scores
	for (size_t i = 0; i < scores.size(); i++) {
		emission_scores[i] /= sum_scores;
	}
}

/*
* conversion of distance scores to emission probabilities using softmin-like function
* given temperature parameter T, emission probability = exp(-score / T) / sum_over_all_alleles(exp(-score / T))

void Entry::convert_scores_to_softmin_probability(unsigned int temperature) {
	assert(scores.size() > 0);
	long double min_score = *std::min_element(scores.begin(), scores.end());
	long double sum_scores = 0.0L;
	emission_scores.resize(scores.size());
	for (size_t i = 0; i < scores.size(); i++) {
		emission_scores[i] = exp(-((long double)scores[i] - min_score) / (long double)temperature);
		sum_scores += emission_scores[i];
	}
	// normalizing the emission scores
	for (size_t i = 0; i < scores.size(); i++) {
		emission_scores[i] /= sum_scores;
	}
}
*/

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

