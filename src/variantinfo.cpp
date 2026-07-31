#include "variantinfo.h"
#include <cassert>

variant_information_t::variant_information_t():
    position(0),
    active_alleles(),
    genotype_likelihoods(),
    allele_references(),
    is_sv(false),
    phasable(false) {}

variant_information_t::variant_information_t(uint32_t pos, uint32_t ploidy, const uint32_t n_alleles, const std::vector<int>& allele_refs, bool sv_flag): 
    position(pos),
    active_alleles(std::vector<bool>(n_alleles, true)),
    genotype_likelihoods(n_alleles, ploidy),
    allele_references(allele_refs),
    is_sv(sv_flag) {
    
    active_gts = std::vector<bool>(genotype_likelihoods.size(), true);
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

std::vector<uint32_t> variant_information_t::get_active_allele_positions() const {
    std::vector<uint32_t> indices;
    for (size_t i = 0; i < active_alleles.size(); i++) {
        if (active_alleles[i]) {
            indices.push_back(i);
        }
    }
    return indices;
}

uint32_t variant_information_t::count_active_gts() const {
    uint32_t count = 0;
    for (size_t i = 0; i < active_gts.size(); i++) {
        if (active_gts[i]) {
            count++;
        }
    }
    return count;
}

std::vector<uint32_t> variant_information_t::get_active_gts_positions() const {
    std::vector<uint32_t> indices;
    for (size_t i = 0; i < active_gts.size(); i++) {
        if (active_gts[i]) {
            indices.push_back(i);
        }
    }
    return indices;
}

void variant_information_t::update_active_alleles(const uint32_t ploidy, const std::vector<uint32_t>& selected_genotype_indices) {
    
    std::fill(active_alleles.begin(), active_alleles.end(), false);
    std::fill(active_gts.begin(), active_gts.end(), false);
    // getting all alleles from selected genotype indices and setting them as true.
    for (auto& selected_genotype_index: selected_genotype_indices) {
        active_gts[selected_genotype_index] = true;
        std::vector<uint32_t> alleles = convert_index_to_alleles(selected_genotype_index, ploidy);
        for (auto allele: alleles) {
            active_alleles[allele] = true;
        }
    }
    if (count_active_alleles() <= 2) {
        assert(count_active_gts() <= 3);
        phasable = true;
    } else {
        phasable = false;
    }
}

void variant_information_t::update_active_gts(const uint32_t ploidy) {
    std::fill(active_gts.begin(), active_gts.end(), false);
    const std::vector<uint32_t> active_alleles = get_active_allele_positions();
    uint32_t allele_i;
    uint32_t allele_j;
    uint32_t gt_index;
    for (int i = 0; i < active_alleles.size(); i++) {
        allele_i = active_alleles[i];
        for (int j = i; j < active_alleles.size(); j++) {
            allele_j = active_alleles[j];
            assert(allele_i <= allele_j);
            gt_index = convert_alleles_to_index(std::vector<uint32_t>{allele_i, allele_j});
            active_gts[gt_index] = true;
        }
    }
}

void variant_information_t::set_as_phasable() {
    assert(count_active_alleles() <= 2);
    assert(count_active_gts() <= 3);
    phasable = true;
}