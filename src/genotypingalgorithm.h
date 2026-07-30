/*
This code controls the alternating phasing-genotyping algorithm.
*/

#ifndef GENOTYPINGALGORITHM_H
#define GENOTYPINGALGORITHM_H

#include "binomial.h"
#include "genotypelikelihoods.h"
#include "variantinfo.h"
#include "readset.h"

#include <cassert>

class GenotypeHMM;
class PhasingDPTable;

class GenotypingAlgorithm {

	public:

		GenotypingAlgorithm(
			ReadSet* read_set, 
			const uint32_t& n_references, 
			const uint32_t& ploidy, 
			const float& temperature, 
			const float& recombrate,
			const float& eff_pop_size,
			const uint32_t& num_sampled_haplotypes,
			const bool remove_reference_path,
			const float& sampling_eff_pop_size,
			const uint32_t& allele_penalty,
			const std::vector<uint32_t>* positions, 
			const std::vector<uint32_t>* n_allele_positions, 
			const std::vector<std::vector<int> >* allele_references, 
			const std::vector<bool>* is_sv_position
		);
		~GenotypingAlgorithm();

		// returns the computed genotype likelihoods for a given position
		std::vector<long double> get_genotype_likelihoods(uint32_t position);

	private:

		friend inline std::ostream& operator<<(std::ostream& out, const GenotypeLikelihoods& g){
			out << "Genotype Likelihoods: ";
			for (size_t i = 0; i < g.as_vector().size(); i++) {
				out << i << ": " << g.as_vector()[i] << ", ";
			}
			return out;
		}

		std::vector<variant_information_t> variant_info_table;
		GenotypeHMM* genotype_hmm;
		PhasingDPTable* phasing_dp_table;

};

#endif
