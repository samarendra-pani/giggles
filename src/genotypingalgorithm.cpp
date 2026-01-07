#include <cassert>

#include "genotypingalgorithm.h"
#include "genotypehmm.h"
#include "phasing/phasingdptable.h"


GenotypingAlgorithm::GenotypingAlgorithm(ReadSet* read_set, const std::vector<float>& recombcost, const uint32_t& n_references, const uint32_t& ploidy, const std::vector<uint32_t>* positions, const std::vector<uint32_t>* n_allele_positions, const std::vector<std::vector<int> >* allele_references, const std::vector<bool>* is_sv_position) {
	
	this->ploidy = ploidy;
	// storing information about each variant position
	variant_info_table = std::vector<variant_information_t>(positions->size());
	//std::vector<GenotypeLikelihoods *> genotype_likelihoods_vector = std::vector<GenotypeLikelihoods *>(positions->size(), nullptr);
	for (size_t i = 0; i < positions->size(); i++) {
		variant_info_table[i] = variant_information_t(positions->at(i), ploidy, n_allele_positions->at(i), allele_references->at(i), is_sv_position->at(i));
		//genotype_likelihoods_vector[i] = &variant_info_table[i].genotype_likelihoods;
	}
	bool is_first_iteration = true;
	read_set->reassignReadIds();
	// some recursive condition to alternate between phasing and genotyping
	while (true) {
		// running the DP table for phasing
		if (ploidy == 2) { 
			/**
			 * PhasingDPTable does the following:
			 * 1. creates the DP table based on the selected reads in the read_set.
			 * 2. creates the optimal haplotypes as "super reads" and the optimal bipartitioning of the reads.
			 * 3. based on the HET positions determined in the super reads, it computes the phasesets.
			 * 4. it haplotags the selected reads based on the optimal bipartitioning.
			 * 5. it haplotags the unselected reads based on the super reads.
			 * If an unselected read does not completely overlap a phaseset, it remains untagged.
			 */
			phasing_dp_table = new PhasingDPTable(read_set, &variant_info_table, is_first_iteration);
			delete phasing_dp_table;
		}
		is_first_iteration = false;
		// running the HMM for genotyping
		exit(0);
		genotype_hmm = new GenotypeHMM(read_set, recombcost, n_references, &variant_info_table);
		// reseting all tags for the next phasing iteration
		read_set->resetTags();
		// break condition to stop alternating
		{
			// some condition
		}
	}
	
}

GenotypingAlgorithm::~GenotypingAlgorithm() {
	delete genotype_hmm;
	
}

std::vector<long double> GenotypingAlgorithm::get_genotype_likelihoods(uint32_t position) {
	// position is the index of the variant position in the VCF file
	return genotype_hmm->get_genotype_likelihoods(position);
}

