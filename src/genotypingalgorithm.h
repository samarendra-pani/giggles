/*
This code controls the alternating phasing-genotyping algorithm.
*/

#ifndef GENOTYPINGALGORITHM_H
#define GENOTYPINGALGORITHM_H

#include "genotypehmm.h"
#include "phasing/phasingdptable.h"
#include "binomial.h"
#include "genotypelikelihoods.h"

#include <cassert>

class GenotypingAlgorithm {

public:
	
	// stores genotyping and phasing information at a given position
	struct variant_information_t {
		unsigned int position; // position in the reference genome
		std::vector<bool> active_alleles; // alleles that were genotyped to be above a certain quality threshold (for use in phasing at further steps)
		std::vector<int> allele_references; // contains the info of which allele came from which reference path (0 - reference, 1... - assemblies/GRCh38)
		GenotypeLikelihoods genotype_likelihoods; // stores genotype likelihoods calculated by the HMM

		variant_information_t() : position(0), active_alleles(), genotype_likelihoods() {}
		variant_information_t(unsigned int pos, const unsigned int n_alleles, const std::vector<int> allele_refs)
			: position(pos), active_alleles(std::vector<bool>(n_alleles, true)), genotype_likelihoods(n_alleles), allele_references(allele_refs) {}

		unsigned int count_active_alleles() const {
			unsigned int count = 0;
			for (size_t i = 0; i < active_alleles.size(); i++) {
				if (active_alleles[i]) {
					count++;
				}
			}
			return count;
		}
	};

	GenotypingAlgorithm(ReadSet* read_set, const std::vector<float>& recombcost, const unsigned int& n_references, const std::vector<unsigned int>* positions, const std::vector<unsigned int>* n_allele_positions, const std::vector<std::vector<int> >* allele_references);
	~GenotypingAlgorithm();

	// returns the computed genotype likelihoods for a given position
	std::vector<long double> get_genotype_likelihoods(unsigned int position);

private:

	friend inline std::ostream& operator<<(std::ostream& out, const GenotypeLikelihoods& g){
		out << "Genotype Likelihoods: ";
		for (int i = 0; i < g.as_vector().size(); i++) {
			out << i << ": " << g.as_vector()[i] << ", ";
		}
		return out;
	}

	friend class GenotypeHMM;
	friend class PhasingDPTable;

	std::vector<variant_information_t> variant_info_table;
	GenotypeHMM* genotype_hmm;
	PhasingDPTable* phasing_dp_table;

};

#endif
