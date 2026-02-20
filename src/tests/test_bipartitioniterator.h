#ifndef TEST_BIPARTITIONITERATOR_H
#define TEST_BIPARTITIONITERATOR_H

/**
 * Unit tests for the Biparition Iterator
 * This builds on top of the Gray Code tests.
 */

#include "../bipartitioniterator.h"
#include "../tests_data.h"
#include <cassert>

void test_clustered_bit_column1(BipartitionIterator* iterator) {
    /**
     * bit index:       | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
     * cluster_id:      | 2 | 3 | 4 | 7 | 8 | 9 |10 |15 |17 |
     * is_clustered:    | T | T | T | T | T | T | T | F | T |
     */
    assert(iterator->is_clustered_bit(0) == true);
    assert(iterator->is_clustered_bit(1) == true);
    assert(iterator->is_clustered_bit(2) == true);
    assert(iterator->is_clustered_bit(3) == true);
    assert(iterator->is_clustered_bit(4) == true);
    assert(iterator->is_clustered_bit(5) == true);
    assert(iterator->is_clustered_bit(6) == true);
    assert(iterator->is_clustered_bit(7) == false);
    assert(iterator->is_clustered_bit(8) == true);

    assert_msg(true, "BipartitionIterator", "Checking for clustered bits passed in column 1.");
}

void test_constrained_bit_column1(BipartitionIterator* iterator) {
    /**
     * bit index:            | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
     * cluster_id:           | 2 | 3 | 4 | 7 | 8 | 9 |10 |15 |17 |
     * has_constrained_bit:  | F | F | F | T | F | T | F | F | F |
     */
    assert(iterator->has_constrained_bit(0) == false);
    assert(iterator->has_constrained_bit(1) == false);
    assert(iterator->has_constrained_bit(2) == false);
    assert(iterator->has_constrained_bit(3) == true);
    assert(iterator->has_constrained_bit(4) == false);
    assert(iterator->has_constrained_bit(5) == true);
    assert(iterator->has_constrained_bit(6) == false);
    assert(iterator->has_constrained_bit(7) == false);
    assert(iterator->has_constrained_bit(8) == false);

    /**
     * bit index:       | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
     * cluster_id:      | 2 | 3 | 4 | 7 | 8 | 9 |10 |15 |17 |
     * constrained_bit: | - | - | - | 0 | - | 2 | - | - | - |
     */
    assert(iterator->get_constrained_bit(3) == 0);
    assert(iterator->get_constrained_bit(5) == 2);

    assert_msg(true, "BipartitionIterator", "Checking for constrained bits in column 1.");
}

void test_advancement_column1(BipartitionIterator* iterator) {
    
    /**
     * Clusters are {2, 3, 4, 7, 8, 9, 10, 15, 17}
     * Free Read Clusters are {2, 3, 4, 8, 10, 15, 17} - These are the ones that will be flipped with Gray Code
     * 
     * Sorting this vector based on number of reads gives {15, 17, 8, 10, 3, 2, 4}
     * 
     * Gray Code ordering happens on this set of cluster where 15 is the LSB and 4 is the MSB -> gives bipartition_index
     * Translating the bits from Gray Code to all the cluster -> gives read_cluster_bit_representation
     */
    
     /**
      * Gray code of size 7
      */
    std::vector<uint32_t> expected_bipartition_indices = {
        0, 1, 3, 2, 6, 7, 5, 4, 12, 13, 15, 14, 10, 11, 9, 8, 
        24, 25, 27, 26, 30, 31, 29, 28, 20, 21, 23, 22, 18, 19, 17, 16, 
        48, 49, 51, 50, 54, 55, 53, 52, 60, 61, 63, 62, 58, 59, 57, 56, 
        40, 41, 43, 42, 46, 47, 45, 44, 36, 37, 39, 38, 34, 35, 33, 32, 
        96, 97, 99, 98, 102, 103, 101, 100, 108, 109, 111, 110, 106, 107, 105, 104, 
        120, 121, 123, 122, 126, 127, 125, 124, 116, 117, 119, 118, 114, 115, 113, 112, 
        80, 81, 83, 82, 86, 87, 85, 84, 92, 93, 95, 94, 90, 91, 89, 88, 
        72, 73, 75, 74, 78, 79, 77, 76, 68, 69, 71, 70, 66, 67, 65, 64
    };
    std::vector<int> expected_bit_changed = { -1,
        7, 8, 7, 4, 7, 8, 7, 6, 7, 8, 7, 4, 7, 8, 7, 1, 
        7, 8, 7, 4, 7, 8, 7, 6, 7, 8, 7, 4, 7, 8, 7, 0, 
        7, 8, 7, 4, 7, 8, 7, 6, 7, 8, 7, 4, 7, 8, 7, 1, 
        7, 8, 7, 4, 7, 8, 7, 6, 7, 8, 7, 4, 7, 8, 7, 2, 
        7, 8, 7, 4, 7, 8, 7, 6, 7, 8, 7, 4, 7, 8, 7, 1, 
        7, 8, 7, 4, 7, 8, 7, 6, 7, 8, 7, 4, 7, 8, 7, 0, 
        7, 8, 7, 4, 7, 8, 7, 6, 7, 8, 7, 4, 7, 8, 7, 1, 
        7, 8, 7, 4, 7, 8, 7, 6, 7, 8, 7, 4, 7, 8, 7
    };
    std::vector<uint32_t> expected_read_cluster_bit_representation = {
        40, 168, 424, 296, 312, 440, 184, 56,
        120, 248, 504, 376, 360, 488, 232, 104,
        106, 234, 490, 362, 378, 506, 250, 122,
        58, 186, 442, 314, 298, 426, 170, 42,
        35, 163, 419, 291, 307, 435, 179, 51,
        115, 243, 499, 371, 355, 483, 227, 99,
        97, 225, 481, 353, 369, 497, 241, 113,
        49, 177, 433, 305, 289, 417, 161, 33,
        5, 133, 389, 261, 277, 405, 149, 21,
        85, 213, 469, 341, 325, 453, 197, 69,
        71, 199, 455, 327, 343, 471, 215, 87,
        23, 151, 407, 279, 263, 391, 135, 7,
        14, 142, 398, 270, 286, 414, 158, 30,
        94, 222, 478, 350, 334, 462, 206, 78,
        76, 204, 460, 332, 348, 476, 220, 92,
        28, 156, 412, 284, 268, 396, 140, 12
    };
    std::vector<std::unordered_map<uint32_t, bool>> expected_changed_reads = {
        {},     // no bits changed for initialization
        {{5, true}},    // bit changed = 7 
        {{7, true}},    // bit changed = 8
        {{5, false}},
        {{2, true}, {6, true}},     // bit changed = 4
        {{5, true}},
        {{7, false}},
        {{5, false}},
        
        {{0, true}, {1, true}},     // bit changed = 6
        
        {{5, true}},
        {{7, true}},
        {{5, false}},
        {{2, false}, {6, false}},
        {{5, true}},
        {{7, false}},
        {{5, false}},

        {{8, true}, {9, true}, {12, true}},     // bit changed = 1

        {{5, true}},    // bit changed = 7 
        {{7, true}},    // bit changed = 8
        {{5, false}},
        {{2, true}, {6, true}},     // bit changed = 4
        {{5, true}},
        {{7, false}},
        {{5, false}},
        
        {{0, false}, {1, false}},     // bit changed = 6
        
        {{5, true}},
        {{7, true}},
        {{5, false}},
        {{2, false}, {6, false}},
        {{5, true}},
        {{7, false}},
        {{5, false}},

        {{14, true}, {3, false}, {4, false}, {13, false}},   // bit changed = 0 (constrained with bit 3)

        {{5, true}},    // bit changed = 7 
        {{7, true}},    // bit changed = 8
        {{5, false}},
        {{2, true}, {6, true}},     // bit changed = 4
        {{5, true}},
        {{7, false}},
        {{5, false}},
        
        {{0, true}, {1, true}},     // bit changed = 6
        
        {{5, true}},
        {{7, true}},
        {{5, false}},
        {{2, false}, {6, false}},
        {{5, true}},
        {{7, false}},
        {{5, false}},

        {{8, false}, {9, false}, {12, false}},     // bit changed = 1

        {{5, true}},    // bit changed = 7 
        {{7, true}},    // bit changed = 8
        {{5, false}},
        {{2, true}, {6, true}},     // bit changed = 4
        {{5, true}},
        {{7, false}},
        {{5, false}},
        
        {{0, false}, {1, false}},     // bit changed = 6
        
        {{5, true}},
        {{7, true}},
        {{5, false}},
        {{2, false}, {6, false}},
        {{5, true}},
        {{7, false}},
        {{5, false}},

        {{15, true}, {17, true}, {10, false}, {11, false}, {16, false}},        // bit changed = 2 (constrained with bit 5)

        {{5, true}},    // bit changed = 7 
        {{7, true}},    // bit changed = 8
        {{5, false}},
        {{2, true}, {6, true}},     // bit changed = 4
        {{5, true}},
        {{7, false}},
        {{5, false}},
        
        {{0, true}, {1, true}},     // bit changed = 6
        
        {{5, true}},
        {{7, true}},
        {{5, false}},
        {{2, false}, {6, false}},
        {{5, true}},
        {{7, false}},
        {{5, false}},

        {{8, true}, {9, true}, {12, true}},     // bit changed = 1

        {{5, true}},    // bit changed = 7 
        {{7, true}},    // bit changed = 8
        {{5, false}},
        {{2, true}, {6, true}},     // bit changed = 4
        {{5, true}},
        {{7, false}},
        {{5, false}},
        
        {{0, false}, {1, false}},     // bit changed = 6
        
        {{5, true}},
        {{7, true}},
        {{5, false}},
        {{2, false}, {6, false}},
        {{5, true}},
        {{7, false}},
        {{5, false}},

        {{14, false}, {3, true}, {4, true}, {13, true}},   // bit changed = 0 (constrained with bit 3)

        {{5, true}},    // bit changed = 7 
        {{7, true}},    // bit changed = 8
        {{5, false}},
        {{2, true}, {6, true}},     // bit changed = 4
        {{5, true}},
        {{7, false}},
        {{5, false}},
        
        {{0, true}, {1, true}},     // bit changed = 6
        
        {{5, true}},
        {{7, true}},
        {{5, false}},
        {{2, false}, {6, false}},
        {{5, true}},
        {{7, false}},
        {{5, false}},

        {{8, false}, {9, false}, {12, false}},     // bit changed = 1

        {{5, true}},    // bit changed = 7 
        {{7, true}},    // bit changed = 8
        {{5, false}},
        {{2, true}, {6, true}},     // bit changed = 4
        {{5, true}},
        {{7, false}},
        {{5, false}},
        
        {{0, false}, {1, false}},     // bit changed = 6
        
        {{5, true}},
        {{7, true}},
        {{5, false}},
        {{2, false}, {6, false}},
        {{5, true}},
        {{7, false}},
        {{5, false}},
    };

    uint32_t count = 0;
    std::unordered_map<uint32_t, bool> changed_reads;
    while (iterator->has_next()) {
        int bit_changed = -1;
        iterator->advance(&bit_changed);
        assert(iterator->get_bipartition_index() == expected_bipartition_indices[count]);
        assert(iterator->get_read_cluster_bit_representation() == expected_read_cluster_bit_representation[count]);
        assert(bit_changed == expected_bit_changed[count]);
        changed_reads.clear();
        iterator->get_changed_reads(bit_changed, changed_reads);
        assert(changed_reads.size() == expected_changed_reads[count].size());
        for (auto const& [read_index, new_bit] : expected_changed_reads[count]) {
            assert(changed_reads.count(read_index) > 0);
            assert(changed_reads[read_index] == new_bit);
        }
        count += 1;
    }
    assert(count == expected_bipartition_indices.size());
    assert_msg(true, "BipartitionIterator", "Advance tests passed in column 1.");
}

void test_clustered_bit_column2(BipartitionIterator* iterator) {
    /**
     * bit index:       | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
     * cluster_id:      | 2 | 3 | 4 | 7 | 8 | 9 |17 |29 |
     * is_clustered:    | T | T | T | T | T | T | T | T |
     */
    assert(iterator->is_clustered_bit(0) == true);
    assert(iterator->is_clustered_bit(1) == true);
    assert(iterator->is_clustered_bit(2) == true);
    assert(iterator->is_clustered_bit(3) == true);
    assert(iterator->is_clustered_bit(4) == true);
    assert(iterator->is_clustered_bit(5) == true);
    assert(iterator->is_clustered_bit(6) == true);
    assert(iterator->is_clustered_bit(7) == true);
    assert_msg(true, "BipartitionIterator", "Checking for clustered bits passed in column 2.");
}

void test_constrained_bit_column2(BipartitionIterator* iterator) {
    /**
     * bit index:              | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
     * cluster_id:             | 2 | 3 | 4 | 7 | 8 | 9 |17 |29 |
     * has_constrained_bit:    | F | F | F | T | F | T | F | T |
     */
    assert(iterator->is_clustered_bit(0) == false);
    assert(iterator->is_clustered_bit(1) == false);
    assert(iterator->is_clustered_bit(2) == false);
    assert(iterator->is_clustered_bit(3) == true);
    assert(iterator->is_clustered_bit(4) == false);
    assert(iterator->is_clustered_bit(5) == true);
    assert(iterator->is_clustered_bit(6) == false);
    assert(iterator->is_clustered_bit(7) == true);

    /**
     * bit index:       | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
     * cluster_id:      | 2 | 3 | 4 | 7 | 8 | 9 |17 |29 |
     * constrained_bit: | - | - | - | 0 | - | 2 | - | 6 |
     */
    assert(iterator->get_constrained_bit(3) == 0);
    assert(iterator->get_constrained_bit(5) == 2);
    assert(iterator->get_constrained_bit(7) == 6);
    assert_msg(true, "BipartitionIterator", "Checking for constrained bits in column 2.");
}


void test_advancement_column2(BipartitionIterator* iterator) {
    
    std::vector<uint32_t> expected_bipartition_indices = {
        0, 1, 3, 2, 6, 7, 5, 4,
        12, 13, 15, 14, 10, 11, 9, 8,
        24, 25, 27, 26, 30, 31, 29, 28,
        20, 21, 23, 22, 18, 19, 17, 16};
    std::vector<int> expected_bit_changed = {-1, 
        0, 1, 0, 4, 0, 1, 0, 6,
        0, 1, 0, 4, 0, 1, 0, 2,
        0, 1, 0, 4, 0, 1, 0, 6,
        0, 1, 0, 4, 0, 1, 0};
    std::vector<uint32_t> expected_read_cluster_bit_representation = {
        168, 161, 163, 170, 186, 179, 177, 184,
        120, 113, 115, 122, 106, 99, 97, 104,
        76, 69, 71, 78, 94, 87, 85, 92,
        156, 149, 151, 158, 142, 135, 133, 140
    };
    std::vector<std::unordered_map<uint32_t, bool>> expected_changed_reads = {
        {{8, true}, {1, false}},      // bit changed = 0
        {{4, true}, {7, true}},       // bit changed = 1
        {{8, false}, {1, true}},
        {{0, true}, {2, true}},     // bit changed = 4
        {{8, true}, {1, false}},
        {{4, false}, {7, false}},
        {{8, false}, {1, true}},
        
        {{3, true}, {12, true}, {13, false}},    // bit changed = 6

        {{8, true}, {1, false}},      // bit changed = 0
        {{4, true}, {7, true}},       // bit changed = 1
        {{8, false}, {1, true}},
        {{0, true}, {2, true}},     // bit changed = 4
        {{8, true}, {1, false}},
        {{4, false}, {7, false}},
        {{8, false}, {1, true}},

        {{9, true}, {11, true}, {5, false}, {6, false}, {10, false}},        // bit changed = 2

        {{8, true}, {1, false}},      // bit changed = 0
        {{4, true}, {7, true}},       // bit changed = 1
        {{8, false}, {1, true}},
        {{0, true}, {2, true}},     // bit changed = 4
        {{8, true}, {1, false}},
        {{4, false}, {7, false}},
        {{8, false}, {1, true}},
        
        {{3, false}, {12, false}, {13, true}},    // bit changed = 6

        {{8, true}, {1, false}},      // bit changed = 0
        {{4, true}, {7, true}},       // bit changed = 1
        {{8, false}, {1, true}},
        {{0, true}, {2, true}},     // bit changed = 4
        {{8, true}, {1, false}},
        {{4, false}, {7, false}},
        {{8, false}, {1, true}},
    };

    uint32_t count = 0;
    std::unordered_map<uint32_t, bool> changed_reads;
    while (iterator->has_next()) {
        int bit_changed = -1;
        iterator->advance(&bit_changed);
        assert(iterator->get_bipartition_index() == expected_bipartition_indices[count]);
        assert(iterator->get_read_cluster_bit_representation() == expected_read_cluster_bit_representation[count]);
        assert(bit_changed == expected_bit_changed[count]);
        changed_reads.clear();
        iterator->get_changed_reads(bit_changed, changed_reads);
        assert(changed_reads.size() == expected_changed_reads[count].size());
        for (auto const& [read_index, new_bit] : expected_changed_reads[count]) {
            assert(changed_reads.count(read_index) > 0);
            assert(changed_reads[read_index] == new_bit);
        }
        count += 1;
    }
    assert(count == expected_bipartition_indices.size());
    assert_msg(true, "BipartitionIterator", "Advance tests passed in column 2.");
}


void test_bipartitioniterator() {
    
    ReadSet* read_set = mock_readset_1();
    std::vector<uint32_t> curr_read_ids = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27};
    std::vector<uint32_t> next_read_ids = {12, 14, 16, 17, 19, 20, 21, 22, 24, 25, 26, 27, 28, 29};
    
    // Testing for first column
    Column* column = new Column(curr_read_ids, next_read_ids, read_set);
    BipartitionIterator* iterator = new BipartitionIterator(column, read_set);

    test_clustered_bit_column1(iterator);
    test_constrained_bit_column1(iterator);
    test_advancement_column1(iterator);

    delete column;
    delete iterator;

    // testing for second column
    Column* column = new Column(next_read_ids, {}, read_set);
    BipartitionIterator* iterator = new BipartitionIterator(column, read_set);

    test_clustered_bit_column2(iterator);
    test_constrained_bit_column2(iterator);
    test_advancement_column2(iterator);

    delete column;
    delete iterator;

    delete read_set;
}


#endif // TEST_BIPARTITIONITERATOR_H