// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#include <cassert>
#include <bits/stdc++.h>
#include <cmath>
using namespace std;

#include "entry.h"

std::vector<long double> Entry::probability_cache;
std::vector<long double> Entry::reciprocal_probability_cache;
uint8_t Entry::k = 100; // less than 256

Entry::Entry(uint32_t r, const std::vector<float>& s) : read_id(r), allele(BLANK) {
	set_scores(s);
}

Entry::Entry(uint32_t r) : read_id(r), scores({}), allele(BLANK) {}

Entry::Entry() : read_id(0), scores({}), allele(BLANK) {}

void Entry::set_read_id(uint32_t r) {
	read_id = r;
}

void Entry::set_scores(const std::vector<float>& s) {
	scores.resize(s.size());
	for (uint32_t i = 0; i < s.size(); i++) {
		float score = s[i];
		uint8_t discretized_score;
		discretized_score = (int)(score*k + 0.5);
		scores[i] = discretized_score;
	}
}

void Entry::set_scores(const std::vector<uint8_t>& s) {
	scores = s;
}

void Entry::set_allele_type(const std::vector<bool>& active_alleles) {

	if (active_alleles.size() > get_max_genotype_alleles()) {
		throw std::runtime_error("Number of alleles greater than max alleles supported.");
	}
	assert(active_alleles.size() == scores.size());
	std::vector<uint8_t> active_scores;
	for (size_t i = 0; i < active_alleles.size(); i++) {
		if (active_alleles[i]) {
			active_scores.push_back(this->scores[i]);
		}
	}
	if (active_scores.size() == 1) { allele = ALLELE1; return; }
	assert(active_scores.size() == 2);
	if (active_scores[0] > active_scores[1]) { allele = ALLELE1; }
	if (active_scores[0] < active_scores[1]) { allele = ALLELE2; }
	if (active_scores[0] == active_scores[1]) { allele = EQUAL_SCORES; }
}

void Entry::set_allele_type(allele_t a) {
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

std::vector<uint8_t> Entry::get_scores() const {
	return scores;
}
		
bool Entry::has_allele_type() const {
	return (allele != BLANK);
}

void Entry::initialize_probability_cache(float temperature) {
    assert (temperature > 0.0f);
    if (!probability_cache.empty()) return;
	
    probability_cache.resize(k+1);
	long double lower_bound = 1e-30L;	// prob(g = 0) = 10^-10
    probability_cache[0] = lower_bound;
	reciprocal_probability_cache.resize(k+1);
    reciprocal_probability_cache[0] = 1.0L/lower_bound;

	long double alpha = ((long double)std::exp(temperature)*lower_bound - 1.0L)/((long double)std::exp(temperature) - 1.0L);
	long double beta = (1.0L - lower_bound)/(std::exp(temperature) - 1.0L);

	for (uint32_t i = 1; i < k; i++) {
		long double g = (1.0L/k) * i;
		probability_cache[i] = (long double)(alpha + (beta * std::exp(g*temperature)));
		reciprocal_probability_cache[i] = 1.0L/probability_cache[i];
	}
	probability_cache[k] = 1.0L;	// prob(g = 1) = 1
	reciprocal_probability_cache[k] = 1.0L;
}

long double Entry::get_emission_score(uint32_t i) const {
	uint8_t score = scores[i];
	return probability_cache[score];
}

long double Entry::get_reciprocal_emission_score(uint32_t i) const {
	uint8_t score = scores[i];
	return reciprocal_probability_cache[score];
}

std::ostream& operator<<(std::ostream& out, const Entry& e) {
	out << "Entry(Read ID: " << e.read_id ;
	out << ", Allele: "<< (int)e.allele;
	out << ")" << std::endl;
	return out;
}

