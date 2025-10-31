#include <cassert>

#include "genotypingalgorithm.h"

GenotypingAlgorithm::GenotypingAlgorithm(ReadSet* read_set, const std::vector<float>& recombcost, const unsigned int& n_references, const std::vector<unsigned int>* positions, const std::vector<unsigned int>* n_allele_positions, const std::vector<std::vector<int> >* allele_references) {
	
	// storing information about each variant position
	variant_info_table = std::vector<variant_information_t>(positions->size());
	//std::vector<GenotypeLikelihoods *> genotype_likelihoods_vector = std::vector<GenotypeLikelihoods *>(positions->size(), nullptr);
	for (size_t i = 0; i < positions->size(); i++) {
		variant_info_table[i] = variant_information_t(positions->at(i), n_allele_positions->at(i), allele_references->at(i));
		//genotype_likelihoods_vector[i] = &variant_info_table[i].genotype_likelihoods;
	}

	// some recursive condition to alternate between phasing and genotyping
	while (true) {
		// running the DP table for phasing
		phasing_dp_table = new PhasingDPTable(read_set, &variant_info_table);
		// running the HMM for genotyping
		genotype_hmm = new GenotypeHMM(read_set, recombcost, n_references, &variant_info_table);
		// break condition to stop alternating
		{
			// some condition
		}
	}
	
}

GenotypingAlgorithm::~GenotypingAlgorithm() {
	delete phasing_dp_table;
	delete genotype_hmm;
	
}

std::vector<long double> GenotypingAlgorithm::get_genotype_likelihoods(unsigned int position) {
	// position is the index of the variant position in the VCF file
	return genotype_hmm->get_genotype_likelihoods(position);
}

