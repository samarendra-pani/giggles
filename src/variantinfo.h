#ifndef VARIANTINFO_H
#define VARIANTINFO_H

#include <vector>
#include <iostream>
#include <numeric>
#include <cstdint>

#include "genotypelikelihoods.h"

struct variant_information_t {

    /**
     * position in the reference genome
     */
    uint32_t position;
    /**
     * alleles that were genotyped to be above a certain quality threshold (for use in phasing at further steps)
     */
    std::vector<bool> active_alleles;
    /**
     * genotypes that have likelihoods above a threshold selected by GenotypeLikelihoods.select_genotypes();
     */
    std::vector<bool> active_gts;
    /**
     * contains the info of which allele came from which reference path (0 - reference, 1... - assemblies/GRCh38)
     */
    std::vector<int> allele_references;
    /**
     * is this variant a structural variant?
     */
    bool is_sv;
    /**
     * Can this variant position be phased?
     * Yes if it has 2 or less active alleles.
     */
    bool phasable; 
    /**
     * stores genotype likelihoods calculated by the HMM
     */
    GenotypeLikelihoods genotype_likelihoods;

    variant_information_t();
    variant_information_t(uint32_t pos, uint32_t ploidy, const uint32_t n_alleles, const std::vector<int>& allele_refs, bool sv_flag);
    
    // get number of alleles defined at this position
    uint32_t get_num_alleles() const;
    
    // count number of active alleles
    uint32_t count_active_alleles() const;

    // count number of active gts
    uint32_t count_active_gts() const;

    // get active position indices
    std::vector<uint32_t> get_active_allele_positions() const;

    // get active gts indices
    std::vector<uint32_t> get_active_gts_positions() const;

    /**
     * Update active alleles based on the genotype likelihoods
     */
    void update_active_alleles(const uint32_t ploidy, const std::vector<uint32_t>& selected_genotype_indices);

    // set active gts based on active alleles (used for active alleles out of HaplotypeSampler)
    void update_active_gts(const uint32_t ploidy);

    // setting the variant position as phasable
    void set_as_phasable();
};

#endif // VARIANTINFO_H