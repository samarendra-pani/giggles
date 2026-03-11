#include "test_genotype.h"

void test_empty_genotype() {
    Genotype g;
    assert_msg(g.is_none(), "Genotype", "Default constructor should be none");
    assert_msg(g.get_ploidy() == 0, "Genotype", "Default constructor should have ploidy 0");
    assert_msg(g.toString() == ".", "Genotype", "Empty genotype string should be '.'");
}

void test_diploid_construction_and_ordering() {
    // Test 0/0
    Genotype g00(std::vector<uint32_t>{0, 0});
    assert_msg(g00.get_ploidy() == 2, "Genotype", "Ploidy should be 2");
    assert_msg(g00.is_homozygous(), "Genotype", "0/0 should be homozygous");
    assert_msg(g00.toString() == "0/0", "Genotype", "String rep should be 0/0");
    assert_msg(g00.get_index() == 0, "Genotype", "Canonical index of 0/0 should be 0");

    // Test 0/1 (Sorting check)
    Genotype g01(std::vector<uint32_t>{0, 1});
    Genotype g10(std::vector<uint32_t>{1, 0});
    
    assert_msg(g01 == g10, "Genotype", "Genotypes {0,1} and {1,0} should be equal (internal sorting)");
    assert_msg(g01.toString() == "0/1", "Genotype", "String rep should be ascending (0/1)");
    assert_msg(!g01.is_homozygous(), "Genotype", "0/1 is not homozygous");
    
    // Note: Based on your code, as_vector returns descending order [Large, Small]
    // because position 0 holds the largest allele.
    std::vector<uint32_t> vec = g01.as_vector();
    assert_msg(vec[0] == 1 && vec[1] == 0, "Genotype", "as_vector() should return descending order [1, 0]");
}

void test_canonical_index_mapping() {
    /** 
     * Mappings for Ploidy 2:
     * 0 -> 0/0
     * 1 -> 0/1
     * 2 -> 1/1
     * 3 -> 0/2
     * 4 -> 1/2
     * 5 -> 2/2
     */
    
    Genotype g_idx1(1, 2); // Index 1, Ploidy 2
    assert_msg(g_idx1.toString() == "0/1", "Genotype", "Index 1 should be 0/1");
    
    Genotype g_idx3(3, 2); // Index 3, Ploidy 2
    assert_msg(g_idx3.toString() == "0/2", "Genotype", "Index 3 should be 0/2");

    Genotype g_idx5(5, 2);
    assert_msg(g_idx5.toString() == "2/2", "Genotype", "Index 5 should be 2/2");
    assert_msg(g_idx5.get_index() == 5, "Genotype", "Round trip index check failed for 5");

    // Test conversion free function directly
    std::vector<uint32_t> alleles = convert_index_to_alleles(4, 2); // Index 4 -> 1/2
    // The conversion function returns unsorted raw alleles usually, but let's check values
    bool has_1 = (alleles[0] == 1 || alleles[1] == 1);
    bool has_2 = (alleles[0] == 2 || alleles[1] == 2);
    assert_msg(has_1 && has_2, "Genotype", "convert_index_to_alleles(4,2) should yield 1 and 2");
}

void test_haploid_chrY() {
    // Testing Ploidy 1
    Genotype g_hap(std::vector<uint32_t>{5});
    assert_msg(g_hap.get_ploidy() == 1, "Genotype", "Ploidy should be 1");
    assert_msg(g_hap.toString() == "5", "Genotype", "String rep should be '5'");
    
    // Check index calculation for haploid
    // For k=1: index = binomial(1 + allele - 1, allele - 1) ... wait, standard formula implies:
    // 0->0, 1->1, 2->2 for haploid.
    assert_msg(g_hap.get_index() == 5, "Genotype", "Haploid index should equal allele number");
}

void test_predicates() {
    Genotype g_bi(std::vector<uint32_t>{0, 1});
    assert_msg(g_bi.is_diploid_and_biallelic(), "Genotype", "0/1 should be diploid biallelic");
    
    Genotype g_multi(std::vector<uint32_t>{0, 2});
    assert_msg(!g_multi.is_diploid_and_biallelic(), "Genotype", "0/2 contains allele > 1");
    
    Genotype g_hom(std::vector<uint32_t>{1, 1});
    assert_msg(g_hom.is_homozygous(), "Genotype", "1/1 is homozygous");
    assert_msg(g_hom.is_diploid_and_biallelic(), "Genotype", "1/1 is biallelic (alleles <= 1)");
}

void test_limits_and_exceptions() {
    // 1. Test Max Ploidy Exceeded
    try {
        Genotype g(std::vector<uint32_t>{0, 1, 2});
        assert_msg(false, "Genotype", "Should have thrown runtime_error for ploidy > 2");
    } catch (const std::runtime_error& e) {
        // Expected
        assert_msg(true, "Genotype", "Should have thrown runtime_error for ploidy > 2");
    }

    // 2. Test Max Allele Exceeded
    try {
        Genotype g(std::vector<uint32_t>{0, 32768});
        assert_msg(false, "Genotype", "Should have thrown runtime_error for allele >= 32768");
    } catch (const std::runtime_error& e) {
        // Expected
        assert_msg(true, "Genotype", "Should have thrown runtime_error for allele >= 32768");
    }
}

void test_comparison_operators() {
    Genotype g1(1, 2); // 0/1
    Genotype g2(1, 2); // 0/1
    Genotype g3(2, 2); // 1/1
    
    assert_msg(g1 == g2, "Genotype", "Equality operator failed");
    assert_msg(g1 != g3, "Genotype", "Inequality operator failed");
    assert_msg(g1 < g3, "Genotype", "Less than operator failed (index 1 < index 2)");
}

void test_genotype() {

    test_empty_genotype();
    test_diploid_construction_and_ordering();
    test_canonical_index_mapping();
    test_haploid_chrY();
    test_predicates();
    test_limits_and_exceptions();
    test_comparison_operators();

}