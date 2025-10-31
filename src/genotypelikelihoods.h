/*
Taken from whatshap (version 2.8)
Original filename: src/phredgenotypelikelihoods.h
*/

#ifndef GENOTYPE_LIKELIHOODS_H
#define GENOTYPE_LIKELIHOODS_H

#include <array>
#include "genotype.h"

class GenotypeLikelihoods {
public:
	GenotypeLikelihoods(const std::vector<long double>& gl, unsigned int num_alleles);
	GenotypeLikelihoods(unsigned int num_alleles);
	GenotypeLikelihoods();

	long double get_by_genotype(Genotype genotype) const;	// get likelihood for given genotype
	void set_by_genotype(Genotype genotype, long double value);	// set likelihood for given genotype

	long double get_by_index(unsigned int index) const;	// get likelihood for given index
	void set_by_index(unsigned int index, long double value);	// set likelihood for given index

	void increment_by_index(unsigned int index, long double value);	// increment likelihood for given index

	std::string toString() const;	// string representation

	unsigned int get_num_alleles() const;	// get number of alleles

	unsigned int size() const;

	const std::vector<long double>& as_vector() const;

	void get_genotypes(std::vector<Genotype>& genotypes) const;

	std::vector<unsigned int> getPhredScores() const;
	unsigned int getPhredScore(Genotype genotype) const;

	void divide_likelihoods_by(long double& val);

private:
	std::vector<long double> gl;
	unsigned int num_alleles;
};


#endif
