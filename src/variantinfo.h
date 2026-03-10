#ifndef VARIANTINFO_H
#define VARIANTINFO_H

#include <vector>
#include <iostream>
#include <numeric>

#include "genotypelikelihoods.h"

struct variant_information_t {
    uint32_t position; // position in the reference genome
    std::vector<bool> active_alleles; // alleles that were genotyped to be above a certain quality threshold (for use in phasing at further steps)
    std::vector<int> allele_references; // contains the info of which allele came from which reference path (0 - reference, 1... - assemblies/GRCh38)
    bool is_sv; // is this variant a structural variant?
    GenotypeLikelihoods genotype_likelihoods; // stores genotype likelihoods calculated by the HMM

    variant_information_t() : position(0), active_alleles(), genotype_likelihoods(), allele_references(), is_sv(false) {}
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

    // update active alleles based on genotype likelihoods
    // This function is based on HaplotypeMapper::select_genotypes()
    void update_active_alleles(uint32_t ploidy) {
        std::vector<std::pair<long double, uint32_t>> indexed_values;
        uint32_t gl_size = genotype_likelihoods.size();
        indexed_values.reserve(gl_size);

        long double sum = 0.0L;
        for (uint32_t i = 0; i < gl_size; i++) {
            long double value = genotype_likelihoods.get_by_index(i);
            indexed_values.push_back({value, i});
            sum += value;
        }

        if (sum == 0.0L) {
            /**
             * no genotype likelihood values available yet.
             * returning all genotypes as possible
             */
            return;
        }

        /**
         * sorting the genotype likelihoods in descending order.
         */
        std::sort(indexed_values.begin(), indexed_values.end(), 
            [](const std::pair<long double, uint32_t>& a, const std::pair<long double, uint32_t>& b) {
                return a.first > b.first; 
            }
        );

        std::vector<uint32_t> selected_genotype_indices;
        long double current_sum = 0.0L;
        const long double THRESHOLD = 0.9L;

        uint32_t count = 0;
        for (const auto& item : indexed_values) {
            current_sum += item.first;
            selected_genotype_indices.push_back(item.second);
            count += 1;
            /**
             * @todo
             * How to handle cases when the genotyping has extreme confidence in a homozygous variant (so there is only one active allele)?
             * On subsequent phasing, we need two active alleles.
             * 
             * Right now, I enqsure that two genotypes are selected which will ensure at least 2 alleles.
             */
            if (current_sum > THRESHOLD && selected_genotype_indices.size() > 1) {
                break;
            }
        }
        // setting all alleles as false
        std::fill(active_alleles.begin(), active_alleles.end(), false);
        // getting all alleles from selected genotype indices and setting them as true.
        for (auto& selected_genotype_index: selected_genotype_indices) {
            std::vector<uint32_t> alleles = convert_index_to_alleles(selected_genotype_index, ploidy);
            for (auto allele: alleles) {
                active_alleles[allele] = true;
            }
        }
    }
};

#endif // VARIANTINFO_H