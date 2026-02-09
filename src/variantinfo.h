#ifndef VARIANTINFO_H
#define VARIANTINFO_H

#include <vector>
#include <iostream>
#include "genotypelikelihoods.h"

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

    // get active position indices
    std::vector<uint32_t> get_active_positions() const {
        std::vector<uint32_t> indices;
        for (size_t i = 0; i < active_alleles.size(); i++) {
            if (active_alleles[i]) {
                indices.push_back(i);
            }
        }
        return indices;
    }

    // set allele to inactive
    void set_allele_inactive(uint32_t allele_index) {
        if (allele_index < active_alleles.size()) {
            active_alleles[allele_index] = false;
        } else {
            throw std::runtime_error("variant_information_t::set_allele_inactive: allele_index out of bounds.");
        }
    }
};

#endif // VARIANTINFO_H