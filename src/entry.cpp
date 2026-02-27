// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#include <cassert>
#include <bits/stdc++.h>
#include <cmath>
using namespace std;

#include "entry.h"

Entry::Entry(uint32_t r, const std::vector<uint32_t>& s) : 
	read_id(r), allele(BLANK), gt(Genotype()) {
		convert_scores_to_probability(s);
	}

Entry::Entry() : read_id(0), emission_scores({}), allele(BLANK), gt(Genotype()) {}

void Entry::set_read_id(uint32_t r) {
	read_id = r;
}

void Entry::set_scores(const std::vector<uint32_t>& s) {
	convert_scores_to_probability(s);
}

void Entry::set_allele_type(const std::vector<bool>& active_alleles) {

	if (active_alleles.size() > get_max_genotype_alleles()) {
		throw std::runtime_error("Number of alleles greater than max alleles supported.");
	}
	assert(active_alleles.size() == emission_scores.size());
	std::vector<long double> active_scores;
	std::vector<uint32_t> indices;
	for (size_t i = 0; i < active_alleles.size(); i++) {
		if (active_alleles[i]) {
			active_scores.push_back(this->emission_scores[i]);
			indices.push_back(i);
		}
	}
	assert(active_scores.size() == 2);
	gt = Genotype(indices);
	if (active_scores[0] > active_scores[1]) { allele = ALLELE1; }
	if (active_scores[0] < active_scores[1]) { allele = ALLELE2; }
	if (active_scores[0] == active_scores[1]) { allele = EQUAL_SCORES; }
}

void Entry::set_allele_type(allele_t a, uint32_t idx1, uint32_t idx2) {
	if (idx1 >= get_max_genotype_alleles() || idx2 >= get_max_genotype_alleles()) {
		throw std::runtime_error("Number of alleles greater than max alleles supported.");
	}
	allele = a;
	gt = Genotype(std::vector<uint32_t>{idx1, idx2});
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

uint32_t Entry::get_allele() const {
	if (gt.is_none() || !has_allele_type()) {
		return (uint32_t)-2;
	}
	assert(gt.get_ploidy() == 2);
	std::vector<uint32_t> alleles = gt.as_vector();
	/**
	 * Genotype class stores the genotypes in descending order.
	 * So active indices 0 and 2 are stored as 2/0.
	 * 
	 * So, to retrieve ALLELE1, we need to access alleles[1] and vice-versa.
	 */
	switch (get_allele_type()) {
		case ALLELE1:
			return alleles[1];
		case ALLELE2:
			return alleles[0];
		case EQUAL_SCORES:
			return (uint32_t)-1;
		default:
			return (uint32_t)-2;
	}
}

bool Entry::has_allele_type() const {
	return (allele != BLANK);
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
			long double logprob = (long double)std::min(scores[i]*2, (uint32_t)60);
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

