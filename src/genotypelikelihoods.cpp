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

GenotypeLikelihoods::GenotypeLikelihoods(const vector<long double>& gl, uint32_t num_alleles, uint32_t ploidy) : gl(gl), num_alleles(num_alleles) {
	uint32_t expected_size = binomial_coefficient(ploidy + num_alleles - 1, ploidy);
	if (expected_size != this->gl.size()) {
		throw runtime_error("Error: wrong number of given genotype likelihoods given.");
	}
}


GenotypeLikelihoods::GenotypeLikelihoods(uint32_t num_alleles, uint32_t ploidy): num_alleles(num_alleles) {
	this->gl = vector<long double>(binomial_coefficient(ploidy + num_alleles - 1, ploidy), 0.0L);
}


GenotypeLikelihoods::GenotypeLikelihoods() {
	this->num_alleles = 0;
	this->gl = vector<long double>();
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
    oss << ")";
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

void GenotypeLikelihoods::reset() {
	std::fill(gl.begin(), gl.end(), 0.0L);
}

void GenotypeLikelihoods::divide_likelihoods_by(long double& val) {
	std::transform(gl.begin(), gl.end(), gl.begin(), std::bind2nd(std::divides<long double>(), val));
}

std::vector<uint32_t> GenotypeLikelihoods::select_genotypes() const {
    std::vector<std::pair<long double, uint32_t>> indexed_values;
    uint32_t gl_size = gl.size();
    indexed_values.reserve(gl_size);

    long double sum = 0.0L;
    for (uint32_t i = 0; i < gl_size; i++) {
        long double value = get_by_index(i);
        indexed_values.push_back({value, i});
        sum += value;
    }

    if (sum == 0.0L) {
        /**
         * no genotype likelihood values available yet.
         * returning all genotypes as possible
         */
        std::vector<uint32_t> result(gl_size);
        std::iota(std::begin(result), std::end(result), 0);
        return result;
    }

    /**
     * sorting the genotype likelihoods in descending order.
     */
    std::sort(indexed_values.begin(), indexed_values.end(), 
        [](const std::pair<long double, uint32_t>& a, const std::pair<long double, uint32_t>& b) {
            return a.first > b.first; 
        }
    );

    std::vector<uint32_t> result;
    long double current_sum = 0.0L;
    const long double THRESHOLD = 0.9L;

    uint32_t count = 0;
    for (const auto& item : indexed_values) {
        current_sum += item.first;
        result.push_back(item.second);
        count += 1;
        if (current_sum >= THRESHOLD) {
            break;
        }
    }
    return result;     
}