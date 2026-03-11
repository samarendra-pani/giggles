#include "test_haplotypemapper.h"

void test_haplotypemapper() {
    
    std::vector<int> allele_refs = {0, 1, 0, -1, 2, 0, 3, -1, 1, 2};
    variant_information_t var_info = variant_information_t(100, 2, 4, allele_refs, false);
    assert(var_info.genotype_likelihoods.size() == 10);
    
    /**
     * Genotype Likelihoods are all 0. So all genotypes are active.
     */
    HaplotypeMapper mapper(var_info.genotype_likelihoods, var_info.allele_references);
    assert_msg(mapper.get_num_states() == 64, "HaplotypeMapper", "Initial number of states.");
    /**
     * testing get_state_index
     */
    std::vector<std::vector<int>> grid = {
        { 0,  1,  2, -1,  3,  4,  5, -1,  6,  7},
        { 8,  9, 10, -1, 11, 12, 13, -1, 14, 15},
        {16, 17, 18, -1, 19, 20, 21, -1, 22, 23},
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {24, 25, 26, -1, 27, 28, 29, -1, 30, 31},
        {32, 33, 34, -1, 35, 36, 37, -1, 38, 39},
        {40, 41, 42, -1, 43, 44, 45, -1, 46, 47},
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {48, 49, 50, -1, 51, 52, 53, -1, 54, 55},
        {56, 57, 58, -1, 59, 60, 61, -1, 62, 63}
    };
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            assert(mapper.get_state_index(i, j) == grid[i][j]);
        }
    }
    assert_msg(true, "HaplotypeMapper", "Inital state index.");

    // Reverse map: index i -> {row, column}
    std::vector<std::pair<u_int32_t, u_int32_t>> reverse_map = {
        {0, 0}, {0, 1}, {0, 2}, {0, 4}, {0, 5}, {0, 6}, {0, 8}, {0, 9},   // Values 0-7
        {1, 0}, {1, 1}, {1, 2}, {1, 4}, {1, 5}, {1, 6}, {1, 8}, {1, 9},   // Values 8-15
        {2, 0}, {2, 1}, {2, 2}, {2, 4}, {2, 5}, {2, 6}, {2, 8}, {2, 9},   // Values 16-23
        {4, 0}, {4, 1}, {4, 2}, {4, 4}, {4, 5}, {4, 6}, {4, 8}, {4, 9},   // Values 24-31
        {5, 0}, {5, 1}, {5, 2}, {5, 4}, {5, 5}, {5, 6}, {5, 8}, {5, 9},   // Values 32-39
        {6, 0}, {6, 1}, {6, 2}, {6, 4}, {6, 5}, {6, 6}, {6, 8}, {6, 9},   // Values 40-47
        {8, 0}, {8, 1}, {8, 2}, {8, 4}, {8, 5}, {8, 6}, {8, 8}, {8, 9},   // Values 48-55
        {9, 0}, {9, 1}, {9, 2}, {9, 4}, {9, 5}, {9, 6}, {9, 8}, {9, 9}    // Values 56-63
    };
    for (uint32_t i = 0; i < mapper.get_num_states(); i++) {
        assert(mapper.get_haplotypes_indices(i) == reverse_map[i]);
    }
    assert_msg(true, "HaplotypeMapper", "Inital reverse map.");
    
    // corresponding genotypes:             0/0    0/1   1/1   0/2   1/2   2/2   0/3   1/3    2/3    3/3
    std::vector<long double> likelihoods = {0.01L, 0.1L, 0.0L, 0.2L, 0.1L, 0.3L, 0.2L, 0.05L, 0.02L, 0.02L};
    // selcted genotypes are 0/1, 0/2, 1/2, 2/2 and 0/3
    var_info.genotype_likelihoods = GenotypeLikelihoods(likelihoods, 4, 2);
    mapper = HaplotypeMapper(var_info.genotype_likelihoods, var_info.allele_references);

    assert_msg(mapper.get_num_states() == 42, "HaplotypeMapper", "Number of states with selected genotypes.");
    /**
     * testing get_state_index
     * std::vector<int> allele_refs = {0, 1, 0, -1, 2, 0, 3, -1, 1, 2};
     */
    grid = {
        {-1,  0, -1, -1,  1, -1,  2, -1,  3,  4},
        { 5, -1,  6, -1,  7,  8, -1, -1, -1,  9},
        {-1, 10, -1, -1, 11, -1, 12, -1, 13, 14},
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {15, 16, 17, -1, 18, 19, -1, -1, 20, 21},
        {-1, 22, -1, -1, 23, -1, 24, -1, 25, 26},
        {27, -1, 28, -1, -1, 29, -1, -1, -1, -1},
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {30, -1, 31, -1, 32, 33, -1, -1, -1, 34},
        {35, 36, 37, -1, 38, 39, -1, -1, 40, 41}
    };
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            assert(mapper.get_state_index(i, j) == grid[i][j]);
        }
    }
    assert_msg(true, "HaplotypeMapper", "State index with selected genotypes.");

    // Reverse map: index i -> {row, column}
    reverse_map = {
        {0, 1}, {0, 4}, {0, 6}, {0, 8}, {0, 9},
        {1, 0}, {1, 2}, {1, 4}, {1, 5}, {1, 9},
        {2, 1}, {2, 4}, {2, 6}, {2, 8}, {2, 9},
        {4, 0}, {4, 1}, {4, 2}, {4, 4}, {4, 5}, {4, 8}, {4, 9},
        {5, 1}, {5, 4}, {5, 6}, {5, 8}, {5, 9},
        {6, 0}, {6, 2}, {6, 5},
        {8, 0}, {8, 2}, {8, 4}, {8, 5}, {8, 9},
        {9, 0}, {9, 1}, {9, 2}, {9, 4}, {9, 5}, {9, 8}, {9, 9}
    };
    for (uint32_t i = 0; i < mapper.get_num_states(); i++) {
        assert(mapper.get_haplotypes_indices(i) == reverse_map[i]);
    }
    assert_msg(true, "HaplotypeMapper", "Reverse map with selected genotypes.");
}