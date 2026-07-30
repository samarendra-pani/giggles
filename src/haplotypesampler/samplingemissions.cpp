/**
 * Original implementation at https://github.com/eblerjana/pangenie/blob/master/src/samplingemissions.cpp (Commit f682fb6)
 */

#include <map>
#include <cmath>
#include <cassert>
#include <algorithm>
#include "samplingemissions.h"

SamplingEmissions::SamplingEmissions(const std::vector<const Entry*>& entries, uint32_t n_alleles) {
	this->allele_penalties = std::vector<unsigned short>(n_alleles);
	this->default_penalty = 100;
	
	std::vector<unsigned short> allele_frequencies(n_alleles, 0);
	float total = (float)entries.size();

	uint8_t best_score; 
	std::vector<uint8_t> scores;
	for (auto entry : entries) {
		best_score = entry->get_max_score();
		scores = entry->get_scores();
		assert(scores.size() == n_alleles);
		for (int allele = 0; allele < n_alleles; allele++) {
			if (scores[allele] == best_score) {
				allele_frequencies[allele]++;
			}
		}
	}
	int allele=0;
	for (auto a: allele_frequencies) {
		float fraction = (float)a/total;
		if (fraction > 0.0) {
			this->allele_penalties[allele] = -10.0 * log10(fraction);
			assert(this->allele_penalties[allele] < this->default_penalty);
		} else {
			// Note: this value is based on the max number of unique kmers (which is 300)
			this->allele_penalties[allele] = this->default_penalty;
		}
		// std::cout << "Allele" << allele << " EP: " << this->allele_penalties[allele] << std::endl;
		allele++;
	}
}

unsigned int SamplingEmissions::get_emission_cost(unsigned short allele_id) const {
	if (allele_id == (unsigned short)-1) {
		return 200;
	}
	return this->allele_penalties[allele_id];
}

void SamplingEmissions::penalize(unsigned short allele_id, unsigned short penalty) {
	this->allele_penalties[allele_id] += penalty;
	if (this->allele_penalties[allele_id] > this->default_penalty) {
		// make sure max penality value is at most default penalty
		this->allele_penalties[allele_id] = this->default_penalty; // + penalty;
	}
	// std::cout << "Penalizing Allele" << allele_id << " New EP: " << this->allele_penalties[allele_id] << std::endl;
}