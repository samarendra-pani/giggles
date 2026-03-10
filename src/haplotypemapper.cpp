#include <numeric>

#include "haplotypemapper.h"

HaplotypeMapper::HaplotypeMapper(const GenotypeLikelihoods& genotype_likelihoods, const std::vector<int>& allele_references) {

    std::vector<uint32_t> selected_indices = select_genotypes(genotype_likelihoods);
    std::vector<bool> is_genotype_selected(genotype_likelihoods.size(), false);
    for (uint32_t idx : selected_indices) {
        if (idx < is_genotype_selected.size()) {
            is_genotype_selected[idx] = true;
        }
    }
    bool all_genotypes_selected = selected_indices.size() == genotype_likelihoods.size();
    uint32_t n_haplotypes = allele_references.size();
    haplotype_to_selected_map.remake(n_haplotypes, n_haplotypes, -1);
    
    std::vector<uint32_t> sorted_alleles;
    sorted_alleles.reserve(2);
    uint32_t count = 0;
    
    for (uint32_t i = 0; i < n_haplotypes; i++) {
        for (uint32_t j = 0; j < n_haplotypes; j++) {
            int allele_i = allele_references.at(i);
            int allele_j = allele_references.at(j);
            if (allele_i == -1 || allele_j == -1) { continue; } // if one of the alleles is invalid or not available, don't consider that pair.
            sorted_alleles.clear();
            if (allele_i < allele_j) {
                sorted_alleles.push_back((uint32_t)allele_i);
                sorted_alleles.push_back((uint32_t)allele_j);
            }
            else {
                sorted_alleles.push_back((uint32_t)allele_j);
                sorted_alleles.push_back((uint32_t)allele_i);
            }
            uint32_t index = convert_alleles_to_index(sorted_alleles);
            /**
             * index is added if
             *  - all genotypes have been selected.
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

std::vector<uint32_t> HaplotypeMapper::select_genotypes(const GenotypeLikelihoods& genotype_likelihoods) const {

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
        std::vector<uint32_t> result(gl_size);
        std::iota(std::begin(result), std::end(result), 0);
        return result;
    }

    /**
     * sorting the genotype likelihoods in descending order.
     */
    std::sort(indexed_values.begin(), indexed_values.end(), 
        [](const std::pair<long double, uint32_t>& a, const std::pair<long double, uint32_t>& b) {
            return a.first > b.first; 
        }
    );

    std::vector<uint32_t> result;
    long double current_sum = 0.0L;
    const long double THRESHOLD = 0.9L;

    uint32_t count = 0;
    for (const auto& item : indexed_values) {
        current_sum += item.first;
        result.push_back(item.second);
        count += 1;
        if (current_sum >= THRESHOLD) {
            break;
        }
    }
    return result;     
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