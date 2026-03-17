#include "variantinfo.h"

variant_information_t::variant_information_t():
    position(0),
    active_alleles(),
    genotype_likelihoods(),
    allele_references(),
    is_sv(false),
    phasable(false) {}

variant_information_t::variant_information_t(uint32_t pos, uint32_t ploidy, const uint32_t n_alleles, const std::vector<int> allele_refs, bool sv_flag): 
    position(pos),
    active_alleles(std::vector<bool>(n_alleles, true)),
    genotype_likelihoods(n_alleles, ploidy),
    allele_references(allele_refs),
    is_sv(sv_flag) {

    if (n_alleles == 2 && !is_sv) {
        phasable = true;
    }
    else {
        phasable = false;
    }
}

uint32_t variant_information_t::get_num_alleles() const {
    return active_alleles.size();
}

uint32_t variant_information_t::count_active_alleles() const {
    uint32_t count = 0;
    for (size_t i = 0; i < active_alleles.size(); i++) {
        if (active_alleles[i]) {
            count++;
        }
    }
    return count;
}

std::vector<uint32_t> variant_information_t::get_active_positions() const {
    std::vector<uint32_t> indices;
    for (size_t i = 0; i < active_alleles.size(); i++) {
        if (active_alleles[i]) {
            indices.push_back(i);
        }
    }
    return indices;
}

void variant_information_t::update_active_alleles(const uint32_t ploidy, const std::vector<uint32_t>& selected_genotype_indices) {
    
    std::fill(active_alleles.begin(), active_alleles.end(), false);
    // getting all alleles from selected genotype indices and setting them as true.
    for (auto& selected_genotype_index: selected_genotype_indices) {
        std::vector<uint32_t> alleles = convert_index_to_alleles(selected_genotype_index, ploidy);
        for (auto allele: alleles) {
            active_alleles[allele] = true;
        }
    }
    if (count_active_alleles() <= 2) {
        phasable = true;
    } else {
        phasable = false;
    }
}