#include "test_variantinfo.h"

void test_empty_variantinfo() {
    variant_information_t var_info = variant_information_t();
    assert(var_info.position == 0);
    assert(var_info.active_alleles.size() == 0);
    assert(var_info.allele_references.size() == 0);
    assert(var_info.genotype_likelihoods.size() == 0);
    assert(var_info.is_sv == false);

    assert_msg(true, "VariantInfo", "Empty VariantInfo object.");
}

void test_active_allele_functions() {
    std::vector<int> allele_refs = {0, 1, 0, 1, 2, 3, 0};
    variant_information_t var_info = variant_information_t(100, 2, 4, allele_refs, false);

    assert(var_info.get_num_alleles() == 4);
    assert(var_info.count_active_alleles() == var_info.active_alleles.size());

    var_info.active_alleles[0] = false;
    var_info.active_alleles[1] = false;

    assert(var_info.count_active_alleles() == 2);
    std::vector<uint32_t> active_positions = var_info.get_active_positions();
    assert(active_positions.size() == 2);
    assert(active_positions[0] == 2);
    assert(active_positions[1] == 3);

    var_info.active_alleles[3] = false;
    assert(var_info.count_active_alleles() == 1);
    active_positions = var_info.get_active_positions();
    assert(active_positions.size() == 1);
    assert(active_positions[0] == 2);
    
    assert_msg(true, "VariantInfo", "Functions related to active alleles.");
}

void test_genotype_likelihoods() {
    std::vector<int> allele_refs = {0, 1, 0, 1, 2, 3, 0};
    variant_information_t var_info = variant_information_t(100, 2, 4, allele_refs, false);
    
    assert(var_info.genotype_likelihoods.size() == 10);
    std::vector<long double> likelihoods = {0.05L, 0.1L, 0.05L, 0.2L, 0.1L, 0.3L, 0.05L, 0.05L, 0.05L, 0.05L};
    var_info.genotype_likelihoods = GenotypeLikelihoods(likelihoods, 4, 2);

    assert(var_info.genotype_likelihoods.get_num_alleles() == 4);

    assert_msg(true, "VariantInfo", "Basic GenotypeLikelihood test.");
}

void test_update_active_allele() {
    std::vector<int> allele_refs = {0, 1, 0, 1, 2, 3, 0};
    variant_information_t var_info = variant_information_t(100, 2, 3, allele_refs, false);
    std::vector<long double> likelihoods = {0.2L, 0.0L, 0.0L, 0.2L, 0.0L, 0.6L};
    std::vector<uint32_t> selected_genotype_indices;
    var_info.genotype_likelihoods = GenotypeLikelihoods(likelihoods, 3, 2);
    selected_genotype_indices = var_info.genotype_likelihoods.select_genotypes();


    // Indices that are selected are 0, 3, and 5 which correspond to 0/0, 0/2 and 2/2
    assert(var_info.count_active_alleles() == 3);
    var_info.update_active_alleles(2, selected_genotype_indices);
    assert(var_info.count_active_alleles() == 2);
    assert(var_info.active_alleles[0] == true);
    assert(var_info.active_alleles[2] == true);
    
    var_info = variant_information_t(100, 2, 3, allele_refs, false);
    likelihoods = {1.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L};
    var_info.genotype_likelihoods = GenotypeLikelihoods(likelihoods, 3, 2);
    selected_genotype_indices = var_info.genotype_likelihoods.select_genotypes();
    selected_genotype_indices = var_info.genotype_likelihoods.select_genotypes();
    // Indices that are selected are 0 which corresponds 0/0. Since at least two genotypes are needed, the next one is 0/1
    assert(var_info.count_active_alleles() == 3);
    var_info.update_active_alleles(2, selected_genotype_indices);
    assert(var_info.count_active_alleles() == 1);
    assert(var_info.active_alleles[0] == true);
    
    var_info = variant_information_t(100, 2, 3, allele_refs, false);
    likelihoods = {0.0L, 0.0L, 0.0L, 0.0L, 1.0L, 0.0L};
    var_info.genotype_likelihoods = GenotypeLikelihoods(likelihoods, 3, 2);
    selected_genotype_indices = var_info.genotype_likelihoods.select_genotypes();
    // Indices that are selected are 4 which corresponds 1/2. Since at least two genotypes are needed, the next one is 0/0
    assert(var_info.count_active_alleles() == 3);
    var_info.update_active_alleles(2, selected_genotype_indices);
    assert(var_info.count_active_alleles() == 2);
    
    assert_msg(true, "VariantInfo", "Updating active alleles based on genotype likelihoods.");

}

void test_variantinfo() {
    test_empty_variantinfo();
    test_active_allele_functions();
    test_genotype_likelihoods();
    test_update_active_allele();
}