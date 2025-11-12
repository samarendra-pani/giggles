/*
This code controls the alternating phasing-genotyping algorithm.
*/

#ifndef GENOTYPINGALGORITHM_H
#define GENOTYPINGALGORITHM_H

#include "binomial.h"
#include "genotypelikelihoods.h"
#include "readset.h"

#include <cassert>

class GenotypeHMM;
class PhasingDPTable;

class GenotypingAlgorithm {

	public:
		
		// stores genotyping and phasing information at a given position
		struct variant_information_t {
			uint32_t position; // position in the reference genome
			std::vector<bool> active_alleles; // alleles that were genotyped to be above a certain quality threshold (for use in phasing at further steps)
			std::vector<int> allele_references; // contains the info of which allele came from which reference path (0 - reference, 1... - assemblies/GRCh38)
			bool is_sv; // is this variant a structural variant?
			GenotypeLikelihoods genotype_likelihoods; // stores genotype likelihoods calculated by the HMM

			variant_information_t() : position(0), active_alleles(), genotype_likelihoods() {}
			variant_information_t(uint32_t pos, uint32_t ploidy, const uint32_t n_alleles, const std::vector<int> allele_refs, bool sv_flag)
				: position(pos), active_alleles(std::vector<bool>(n_alleles, true)), genotype_likelihoods(n_alleles, ploidy), allele_references(allele_refs), is_sv(sv_flag) {}

			// get number of alleles defined at this position
			uint32_t get_num_alleles() const {
				return active_alleles.size();
			}
			// count number of active alleles
			uint32_t count_active_alleles() const {
				uint32_t count = 0;
				for (size_t i = 0; i < active_alleles.size(); i++) {
					if (active_alleles[i]) {
						count++;
					}
				}
				return count;
			}
		};

		GenotypingAlgorithm(ReadSet* read_set, const std::vector<float>& recombcost, const uint32_t& n_references, const uint32_t& ploidy, const std::vector<uint32_t>* positions, const std::vector<uint32_t>* n_allele_positions, const std::vector<std::vector<int> >* allele_references, const std::vector<bool>* is_sv_position);
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

		uint32_t ploidy;
		std::vector<variant_information_t> variant_info_table;
		GenotypeHMM* genotype_hmm;
		PhasingDPTable* phasing_dp_table;

};

#endif
