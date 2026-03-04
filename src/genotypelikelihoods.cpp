/*
Taken from whatshap (version 2.8)
Original filename: src/phredgenotypelikelihoods.cpp
*/

#include <sstream>
#include <cassert>
#include <cmath>

#include "genotypelikelihoods.h"
#include "binomial.h"

using namespace std;

GenotypeLikelihoods::GenotypeLikelihoods(const vector<long double>& gl, uint32_t num_alleles, uint32_t ploidy) : gl(gl), ploidy(ploidy), num_alleles(num_alleles) {
	uint32_t expected_size = binomial_coefficient(ploidy + num_alleles - 1, ploidy);
	if (expected_size != this->gl.size()) {
		throw runtime_error("Error: wrong number of given genotype likelihoods given.");
	}
}


GenotypeLikelihoods::GenotypeLikelihoods(uint32_t num_alleles, uint32_t ploidy): ploidy(ploidy), num_alleles(num_alleles) {
	this->gl = vector<long double>(binomial_coefficient(ploidy + num_alleles - 1, ploidy), 0.0L);
}


GenotypeLikelihoods::GenotypeLikelihoods() {
	this->num_alleles = 0;
	this->gl = vector<long double>();
}

long double GenotypeLikelihoods::get_by_genotype(Genotype genotype) const {
	uint32_t index = genotype.get_index();
	assert(index < this->gl.size());
	return this->gl[index];
}

void GenotypeLikelihoods::set_by_genotype(Genotype genotype, long double value) {
	uint32_t index = genotype.get_index();
	assert(index < this->gl.size());
	this->gl[index] = value;
}

long double GenotypeLikelihoods::get_by_index(uint32_t index) const {
	assert(index < this->gl.size());
	return this->gl[index];
}

void GenotypeLikelihoods::set_by_index(uint32_t index, long double value) {
	assert(index < this->gl.size());
	this->gl[index] = value;
}

void GenotypeLikelihoods::increment_by_index(uint32_t index, long double value) {
	assert(index < this->gl.size());
	this->gl[index] += value;
}

std::string GenotypeLikelihoods::toString() const {
	ostringstream oss;
	oss << "GenotypeLikelihoods(";
	for (size_t i = 0; i < this->gl.size(); ++i) {
		if (i > 0) oss << ",";
		oss << gl[i];
	}
	return oss.str();
}


uint32_t GenotypeLikelihoods::get_num_alleles() const {
	return this->num_alleles;
}

uint32_t GenotypeLikelihoods::size() const {
	return this->gl.size();
}

const vector<long double>& GenotypeLikelihoods::as_vector() const {
	return this->gl;
}

void GenotypeLikelihoods::get_genotypes(vector<Genotype>& genotypes) const {
	for (uint32_t i = 0; i < this->size(); ++i) {
		genotypes.push_back(Genotype(i, ploidy));
	}
} 


// TODO: need to test this phred score calculation
std::vector<uint32_t> GenotypeLikelihoods::getPhredScores() const {
	long double max = 0.0;
	for (int i=0; i<gl.size(); ++i) {
		if (gl[i] > max) max = gl[i];
	}
	if (max == 0.0) {
		return std::vector<uint32_t>(gl.size(), 0);
	}
	std::vector<uint32_t> phred_scores;
	for (int i=0; i<gl.size(); ++i) {
		long double prob = gl[i] / max;
		uint32_t phred = (uint32_t)(-10.0L * log10(prob));
		phred_scores.push_back(phred);
	}
	return phred_scores;
}

// TODO: need to test this phred score calculation
uint32_t GenotypeLikelihoods::getPhredScore(Genotype genotype) const {
	uint32_t index = genotype.get_index();
	assert(index < this->gl.size());
	long double max = 0.0;
	for (int i=0; i<gl.size(); ++i) {
		if (gl[i] > max) max = gl[i];
	}
	if (max == 0.0) {
		return 0;
	}
	long double prob = gl[index] / max;
	uint32_t phred = (uint32_t)(-10.0L * log10(prob));
	return phred;
}

void GenotypeLikelihoods::divide_likelihoods_by(long double& val) {
	std::transform(gl.begin(), gl.end(), gl.begin(), std::bind2nd(std::divides<long double>(), val));
}