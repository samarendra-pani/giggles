/*
Taken from whatshap (version 2.8)
Original filename: src/phredgenotypelikelihoods.h
*/

#ifndef GENOTYPE_LIKELIHOODS_H
#define GENOTYPE_LIKELIHOODS_H

#include <array>
#include <numeric>
#include "genotype.h"

class GenotypeLikelihoods {
public:
	GenotypeLikelihoods(const std::vector<long double>& gl, uint32_t num_alleles, uint32_t ploidy);
	GenotypeLikelihoods(uint32_t num_alleles, uint32_t ploidy);
	GenotypeLikelihoods();

	long double get_by_index(uint32_t index) const;	// get likelihood for given index
	void set_by_index(uint32_t index, long double value);	// set likelihood for given index

	void increment_by_index(uint32_t index, long double value);	// increment likelihood for given index

	std::string toString() const;	// string representation

	uint32_t get_num_alleles() const;	// get number of alleles

	uint32_t size() const;

	void reset();

	const std::vector<long double>& as_vector() const;

	std::vector<uint32_t> getPhredScores() const;
	uint32_t getPhredScore(Genotype genotype) const;

	void divide_likelihoods_by(long double& val);

	/**
     * Selects the genotypes which are "likely" and returns a vector of their canonical index.
     * The likely genotypes are the m genotypes with the highest likelihoods scores such that sum of their scores >= 0.9
     */
	std::vector<uint32_t> select_genotypes() const;

private:
	std::vector<long double> gl;
	uint32_t num_alleles;
};


#endif
