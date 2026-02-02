// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#include <cassert>
#include <bits/stdc++.h>
#include <cmath>
using namespace std;

#include "entry.h"

Entry::Entry(uint32_t r, const std::vector<uint32_t>& s) : 
	read_id(r), allele(BLANK) {
		convert_scores_to_probability(s);
	}

Entry::Entry() : read_id(0), emission_scores({}), allele(BLANK) {}

void Entry::set_read_id(uint32_t r) {
	read_id = r;
}

void Entry::set_scores(const std::vector<uint32_t>& s) {
	convert_scores_to_probability(s);
}

void Entry::set_allele_type(const std::vector<bool>& active_alleles) {
	assert(active_alleles.size() == emission_scores.size());
	std::vector<uint32_t> active_scores;
	std::vector<uint32_t> indices;
	for (size_t i = 0; i < active_alleles.size(); i++) {
		if (active_alleles[i]) {
			active_scores.push_back(this->emission_scores[i]);
			indices.push_back(i);
		}
	}
	assert(active_scores.size() == 2);
	allele1_idx = indices[0];
	allele2_idx = indices[1];
	if (active_scores[0] > active_scores[1]) { allele = ALLELE1; }
	if (active_scores[0] < active_scores[1]) { allele = ALLELE2; }
	if (active_scores[0] == active_scores[1]) { allele = EQUAL_SCORES; }
}

void Entry::set_allele_type(allele_t a, uint32_t idx1, uint32_t idx2) {
	allele = a;
}

uint32_t Entry::get_read_id() const {
	return read_id;
}

// Currently hardcoded to 30 since Whatshap uses fixed quality scores for phasing.
uint32_t Entry::get_phred_score() const {
	return 30;
}

Entry::allele_t Entry::get_allele_type() const {
	return allele;
}

bool Entry::has_allele_type() const {
	return (allele == BLANK);
}

/*
* conversion of distance scores to emission probabilities using error probability 0.0001
* emission probability = 10^(-max(score*log10(0.0001), 1e-10))
* modelling the probability of observing a read given the true allele and error rate
*/
void Entry::convert_scores_to_probability(const std::vector<uint32_t>& scores) {
	if (scores.size() > 0) {
		long double sum_scores = 0.0L;
		emission_scores.resize(scores.size());
		for (size_t i = 0; i < scores.size(); i++) {
			long double logprob = (long double)std::min(scores[i]*log10(0.0001), 1e-10);
			emission_scores[i] = pow(10.0L, -logprob);
			sum_scores += emission_scores[i];
		}
		// normalizing the emission scores
		for (size_t i = 0; i < scores.size(); i++) {
			emission_scores[i] /= sum_scores;
		}
	}
}

/*
* conversion of distance scores to emission probabilities using softmin-like function
* given temperature parameter T, emission probability = exp(-score / T) / sum_over_all_alleles(exp(-score / T))
*/

/*
void Entry::convert_scores_to_softmin_probability(uint32_t temperature) {
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


std::vector<long double> Entry::get_emission_scores() const {
	return emission_scores;
}

void Entry::set_emission_scores(const std::vector<long double>& scores) {
	emission_scores = scores;
}


std::ostream& operator<<(std::ostream& out, const Entry& e) {
	out << "Entry(" << e.read_id ;
	out << ","<< e.allele << ",(";
	for (auto i : e.emission_scores) {
		out << i << ",";
	}
	out << ")" << std::endl;
	return out;
}

