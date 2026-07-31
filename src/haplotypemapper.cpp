#include <numeric>

#include "haplotypemapper.h"

HaplotypeMapper::HaplotypeMapper(const variant_information_t& variant_info) {

    const std::vector<uint32_t>& selected_indices = variant_info.get_active_gts_positions();
    const std::vector<bool>& is_genotype_selected = variant_info.active_gts;
    bool all_genotypes_selected = selected_indices.size() == is_genotype_selected.size();
    const std::vector<int>& allele_references = variant_info.allele_references;
    uint32_t n_haplotypes = allele_references.size();
    haplotype_to_selected_map.remake(n_haplotypes, n_haplotypes, -1);
    
    std::vector<uint32_t> sorted_alleles(2);
    uint32_t count = 0;
    int allele_i;
    int allele_j;
    for (uint32_t i = 0; i < n_haplotypes; i++) {
        for (uint32_t j = 0; j < n_haplotypes; j++) {
            allele_i = allele_references[i];
            allele_j = allele_references[j];
            if (allele_i == -1 || allele_j == -1) { continue; } // if one of the alleles is invalid or not available, don't consider that pair.
            if (allele_i < allele_j) {
                sorted_alleles[0] = (uint32_t)allele_i;
                sorted_alleles[1] = (uint32_t)allele_j;
            }
            else {
                sorted_alleles[0] = (uint32_t)allele_j;
                sorted_alleles[1] = (uint32_t)allele_i;
            }
            uint32_t index = convert_alleles_to_index(sorted_alleles);
            /**
             * index is added if
             *  - all genotypes have been selected.
             *    OR
             *  - if this was one of the selected genotypes.
             */
            bool allowed = all_genotypes_selected || (index < is_genotype_selected.size() && is_genotype_selected[index]);
            if (allowed) {
                reverse_map.push_back({i, j});
                haplotype_to_selected_map.set(i, j, count);
                count++;
            }
        }
    }
}

uint32_t HaplotypeMapper::get_num_states() const {
    return reverse_map.size();
}

int HaplotypeMapper::get_state_index(uint32_t i, uint32_t j) const {
    assert(i < haplotype_to_selected_map.get_size0());
    assert(j < haplotype_to_selected_map.get_size1());
    return haplotype_to_selected_map.at(i, j);
}

std::pair<uint32_t, uint32_t> HaplotypeMapper::get_haplotypes_indices(uint32_t i) const {
    assert(i < reverse_map.size());
    return reverse_map.at(i);
}