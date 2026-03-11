#include "test_column.h"

void test_read_cluster_ids_column1(Column* column) {
    const std::vector<uint32_t>* read_cluster_ids = column->get_read_cluster_ids();
    std::vector<uint32_t> expected_clusters = {2, 3, 4, 7, 8, 9, 10, 15, 17};
    assert(read_cluster_ids->size() == expected_clusters.size());
    for (size_t i = 0; i < expected_clusters.size(); i++) {
        assert(read_cluster_ids->at(i) == expected_clusters[i]);
    }
    assert_msg(true, "Column", "Read cluster IDs test passed for column 1.");
}

void test_num_bipartitions_column1(Column* column) {
    uint32_t num_bipartitions = column->get_num_bipartition();
    uint32_t expected_num_bipartitions = 1 << 7; // (2, 7) and (4, 9) are constrained
    assert_msg(num_bipartitions == expected_num_bipartitions, "Column", "Number of bipartition from column 1.");
}

void test_constrained_position_map_column1(Column* column) {
    const std::unordered_map<uint32_t, uint32_t>* map = column->get_constrained_position_map();

    /**
     * Clusters are {2, 3, 4, 7, 8, 9, 10, 15, 17}
     *  Constrained Pair (2, 7) -> Position of 7 is 3 and Position of 2 is 0. So map is 0 -> 3.
     *  Constrained Pair (4, 9) -> Position of 9 is 5 and Position of 4 is 2. So map is 2 -> 5.
     */

    assert_msg(map->size() == 2, "Column", "Size of column 1 constrained position map.");
    assert(map->at(0) == 3);
    assert(map->at(2) == 5);
    assert_msg(true, "Column", "Content of column 1 constrained position map.");
}

void test_sorted_free_read_cluster_positions_column1(Column* column) {
    const std::vector<uint32_t>* positions = column->get_sorted_free_read_cluster_positions();
    
    /**
     * Clusters are {2, 3, 4, 7, 8, 9, 10, 15, 17}
     * Free Read Clusters are {2, 3, 4, 8, 10, 15, 17}
     * Number of Reads associated with these clusters are {4, 3, 5, 2, 2, 1, 1}
     * 
     * If we sort out (in ascending order) free read clusters based on how many reads they have we get {15, 17, 8, 10, 3, 2, 4}
     * Their positions in the original read clusters vector are {7, 8, 4, 6, 1, 0, 2}
     */
    
    std::vector<uint32_t> expected = {7, 8, 4, 6, 1, 0, 2};
    assert_msg(positions->size() == expected.size(), "Column", "Size of sorted free read cluster positions of column 1.");
    for (uint32_t i = 0; i < positions->size(); i++) {
        assert(positions->at(i) == expected[i]);
    }
    assert_msg(true, "Column", "Content of sorted free read cluster positions of column 1.");

}

void test_cluster_id_to_read_index_map_column1(Column* column) {
    const std::unordered_map<uint32_t, std::vector<uint32_t>>* cluster_map = column->get_cluster_id_to_read_index_map();
    
    /**
     * Expected mapping:
     *  Cluster 10 -> read indices {0, 1}
     *  Cluster 8 -> read indices {2, 6}
     *  Cluster 7 -> read indices {3, 4, 13}
     *  Cluster 15 -> read indices {5}
     *  Cluster 17 -> read indices {7}
     *  Cluster 3 -> read indices {8, 9, 12}
     *  Cluster 9 -> read indices {10, 11, 16}
     *  Cluster 2 -> read indices {14}
     *  Cluster 4 -> read indices {15, 17}
     */

    assert(cluster_map->size() == 9);
    
    assert(cluster_map->at(10).size() == 2);
    assert(cluster_map->at(10)[0] == 0);
    assert(cluster_map->at(10)[1] == 1);

    assert(cluster_map->at(8).size() == 2);
    assert(cluster_map->at(8)[0] == 2);
    assert(cluster_map->at(8)[1] == 6);
    
    assert(cluster_map->at(7).size() == 3);
    assert(cluster_map->at(7)[0] == 3);
    assert(cluster_map->at(7)[1] == 4);
    assert(cluster_map->at(7)[2] == 13);

    assert(cluster_map->at(15).size() == 1);
    assert(cluster_map->at(15)[0] == 5);

    assert(cluster_map->at(17).size() == 1);
    assert(cluster_map->at(17)[0] == 7);
    
    assert(cluster_map->at(3).size() == 3);
    assert(cluster_map->at(3)[0] == 8);
    assert(cluster_map->at(3)[1] == 9);
    assert(cluster_map->at(3)[2] == 12);

    assert(cluster_map->at(9).size() == 3);
    assert(cluster_map->at(9)[0] == 10);
    assert(cluster_map->at(9)[1] == 11);
    assert(cluster_map->at(9)[2] == 16);

    assert(cluster_map->at(2).size() == 1);
    assert(cluster_map->at(2)[0] == 14);

    assert(cluster_map->at(4).size() == 2);
    assert(cluster_map->at(4)[0] == 15);
    assert(cluster_map->at(4)[1] == 17);

    assert_msg(true, "Column", "Cluster ID to Read IDs map test passed for column 1.");
}

void test_bipartition_compatibility_column1(Column* column) {
    
    /**
     * The input to get_backward_compatible_bipartitions() is the read_cluster_bit_representation of the
     *  read clusters at the next position.
     * Clusters at the next position are:                           | 2 | 3 | 4 | 7 | 8 | 9 | 17 | 29 |
     * Their bit positions in read_cluster_bit_representation are:  | 0 | 1 | 2 | 3 | 4 | 5 | 6  | 7  |
     * Their constrained clusters are:                              | 7 | - | 9 | 2 | - | 4 | 29 | 17 |
     * 
     * So bit pairs (0, 3), (2, 5), and (6, 7) have to be opposites (to denote correct read_cluster_bit_representation).
     * For testing purposes, we will not follow to see if the next_read_cluster_masks are behaving properly. 
     */
    std::vector<uint32_t> compatible_bipartitions; // this will store the compatible bipartitions for all constrained clusters set to 0.
    std::vector<uint32_t> compatible_bipartitions_1; // this will store the compatible bipartitions for all constrained clusters except cluster 4 set to 0.
    std::vector<uint32_t> compatible_bipartitions_2; // this will store the compatible bipartitions for all constrained clusters except cluster 9 set to 0.
    std::vector<uint32_t> compatible_bipartitions_3; // this will store the compatible bipartitions for all constrained clusters except cluster 29 set to 0.
    std::vector<uint32_t> compatible_bipartitions_4; // this will store the compatible bipartitions for all constrained clusters set to 1.
    

    
    /**
     * Setting read_cluster_bit_representation to 0 gives back the cached bipartitions
     */
    column->get_backward_compatible_bipartitions(0, compatible_bipartitions);
    assert(compatible_bipartitions.size() == 4);
    assert(compatible_bipartitions[0] == 0);
    assert(compatible_bipartitions[2] == 1);
    assert(compatible_bipartitions[1] == 8);
    assert(compatible_bipartitions[3] == 9);

    /**
     * The above result should remain the same if we flipped the bits for 7, 9 and 27 (bits 3, 5, an 7)
     * So we flip each bit (and one case where all are flipped).
     * We are not considering all permutations (too many to check)
     * 
     * Clusters at the next position are:                           | 2 | 3 | 4 | 7 | 8 | 9 | 17 | 29 |
     * Their bit positions in read_cluster_bit_representation are:  | 0 | 1 | 2 | 3 | 4 | 5 | 6  | 7  |
     * Bit values:                                                  | 0 | 0 | 0 |0/1| 0 |0/1| 0  |0/1 |
     */
    column->get_backward_compatible_bipartitions(8, compatible_bipartitions_1);
    column->get_backward_compatible_bipartitions(32, compatible_bipartitions_2);
    column->get_backward_compatible_bipartitions(128, compatible_bipartitions_3);
    column->get_backward_compatible_bipartitions(168, compatible_bipartitions_4);
    assert(compatible_bipartitions_1 == compatible_bipartitions_2);
    assert(compatible_bipartitions_3 == compatible_bipartitions_4);
    assert(compatible_bipartitions_1 == compatible_bipartitions_3);
    assert(compatible_bipartitions_1 == compatible_bipartitions);

    /** 
     * Clusters at the next position are:                           | 2 | 3 | 4 | 7 | 8 | 9 | 17 | 29 |
     * Their bit positions in read_cluster_bit_representation are:  | 0 | 1 | 2 | 3 | 4 | 5 | 6  | 7  |
     * Bit values:                                                  | 1 | 0 | 0 |0/1| 0 |0/1| 0  |0/1 |
     * 
     * Mask for read cluster 2 is 1 << 5 = 32
     */
    column->get_backward_compatible_bipartitions(1, compatible_bipartitions);
    
    assert(compatible_bipartitions[0] == 32);
    assert(compatible_bipartitions[2] == 33);
    assert(compatible_bipartitions[1] == 40);
    assert(compatible_bipartitions[3] == 41);
    
    column->get_backward_compatible_bipartitions(9, compatible_bipartitions_1);
    column->get_backward_compatible_bipartitions(33, compatible_bipartitions_2);
    column->get_backward_compatible_bipartitions(129, compatible_bipartitions_3);
    column->get_backward_compatible_bipartitions(169, compatible_bipartitions_4);
    assert(compatible_bipartitions_1 == compatible_bipartitions_2);
    assert(compatible_bipartitions_3 == compatible_bipartitions_4);
    assert(compatible_bipartitions_1 == compatible_bipartitions_3);
    assert(compatible_bipartitions_1 == compatible_bipartitions);

    /** 
     * Clusters at the next position are:                           | 2 | 3 | 4 | 7 | 8 | 9 | 17 | 29 |
     * Their bit positions in read_cluster_bit_representation are:  | 0 | 1 | 2 | 3 | 4 | 5 | 6  | 7  |
     * Bit values:                                                  | 0 | 1 | 0 |0/1| 0 |0/1| 0  |0/1 |
     * 
     * Mask for read cluster 3 is 1 << 4 = 16
     */
    column->get_backward_compatible_bipartitions(2, compatible_bipartitions);
    
    assert(compatible_bipartitions[0] == 16);
    assert(compatible_bipartitions[2] == 17);
    assert(compatible_bipartitions[1] == 24);
    assert(compatible_bipartitions[3] == 25);

    column->get_backward_compatible_bipartitions(10, compatible_bipartitions_1);
    column->get_backward_compatible_bipartitions(34, compatible_bipartitions_2);
    column->get_backward_compatible_bipartitions(130, compatible_bipartitions_3);
    column->get_backward_compatible_bipartitions(170, compatible_bipartitions_4);
    assert(compatible_bipartitions_1 == compatible_bipartitions_2);
    assert(compatible_bipartitions_3 == compatible_bipartitions_4);
    assert(compatible_bipartitions_1 == compatible_bipartitions_3);
    assert(compatible_bipartitions_1 == compatible_bipartitions);
    
    /** 
     * Clusters at the next position are:                           | 2 | 3 | 4 | 7 | 8 | 9 | 17 | 29 |
     * Their bit positions in read_cluster_bit_representation are:  | 0 | 1 | 2 | 3 | 4 | 5 | 6  | 7  |
     * Bit values:                                                  | 0 | 0 | 1 |0/1| 0 |0/1| 0  |0/1 |
     * 
     * Mask for read cluster 4 is 1 << 6 = 64
     */
    column->get_backward_compatible_bipartitions(4, compatible_bipartitions);

    assert(compatible_bipartitions[0] == 64);
    assert(compatible_bipartitions[2] == 65);
    assert(compatible_bipartitions[1] == 72);
    assert(compatible_bipartitions[3] == 73);

    column->get_backward_compatible_bipartitions(12, compatible_bipartitions_1);
    column->get_backward_compatible_bipartitions(36, compatible_bipartitions_2);
    column->get_backward_compatible_bipartitions(132, compatible_bipartitions_3);
    column->get_backward_compatible_bipartitions(172, compatible_bipartitions_4);
    assert(compatible_bipartitions_1 == compatible_bipartitions_2);
    assert(compatible_bipartitions_3 == compatible_bipartitions_4);
    assert(compatible_bipartitions_1 == compatible_bipartitions_3);
    assert(compatible_bipartitions_1 == compatible_bipartitions);

    /** 
     * Clusters at the next position are:                           | 2 | 3 | 4 | 7 | 8 | 9 | 17 | 29 |
     * Their bit positions in read_cluster_bit_representation are:  | 0 | 1 | 2 | 3 | 4 | 5 | 6  | 7  |
     * Bit values:                                                  | 0 | 0 | 0 |0/1| 1 |0/1| 0  |0/1 |
     * 
     * Mask for read cluster 8 is 1 << 2 = 4
     */
    column->get_backward_compatible_bipartitions(16, compatible_bipartitions);
    
    assert(compatible_bipartitions[0] == 4);
    assert(compatible_bipartitions[2] == 5);
    assert(compatible_bipartitions[1] == 12);
    assert(compatible_bipartitions[3] == 13);

    column->get_backward_compatible_bipartitions(24, compatible_bipartitions_1);
    column->get_backward_compatible_bipartitions(48, compatible_bipartitions_2);
    column->get_backward_compatible_bipartitions(144, compatible_bipartitions_3);
    column->get_backward_compatible_bipartitions(184, compatible_bipartitions_4);
    assert(compatible_bipartitions_1 == compatible_bipartitions_2);
    assert(compatible_bipartitions_3 == compatible_bipartitions_4);
    assert(compatible_bipartitions_1 == compatible_bipartitions_3);
    assert(compatible_bipartitions_1 == compatible_bipartitions);

    /** 
     * Clusters at the next position are:                           | 2 | 3 | 4 | 7 | 8 | 9 | 17 | 29 |
     * Their bit positions in read_cluster_bit_representation are:  | 0 | 1 | 2 | 3 | 4 | 5 | 6  | 7  |
     * Bit values:                                                  | 0 | 0 | 0 |0/1| 0 |0/1| 1  |0/1 |
     * 
     * Mask for read cluster 17 is 1 << 1 = 2
     */
    column->get_backward_compatible_bipartitions(64, compatible_bipartitions);

    assert(compatible_bipartitions[0] == 2);
    assert(compatible_bipartitions[2] == 3);
    assert(compatible_bipartitions[1] == 10);
    assert(compatible_bipartitions[3] == 11);

    column->get_backward_compatible_bipartitions(72, compatible_bipartitions_1);
    column->get_backward_compatible_bipartitions(96, compatible_bipartitions_2);
    column->get_backward_compatible_bipartitions(192, compatible_bipartitions_3);
    column->get_backward_compatible_bipartitions(232, compatible_bipartitions_4);
    assert(compatible_bipartitions_1 == compatible_bipartitions_2);
    assert(compatible_bipartitions_3 == compatible_bipartitions_4);
    assert(compatible_bipartitions_1 == compatible_bipartitions_3);
    assert(compatible_bipartitions_1 == compatible_bipartitions);

    assert_msg(true, "Column", "Backward compatibility tests passed.");
}  

void test_read_cluster_ids_column2(Column* column) {
    const std::vector<uint32_t>* read_cluster_ids = column->get_read_cluster_ids();
    std::vector<uint32_t> expected_clusters = {2, 3, 4, 7, 8, 9, 17, 29};
    assert(read_cluster_ids->size() == expected_clusters.size());
    for (size_t i = 0; i < expected_clusters.size(); i++) {
        assert(read_cluster_ids->at(i) == expected_clusters[i]);
    }
    assert_msg(true, "Column", "Read cluster IDs test passed for column 2.");
}

void test_num_bipartitions_column2(Column* column) {
    uint32_t num_bipartitions = column->get_num_bipartition();
    uint32_t expected_num_bipartitions = 1 << 5; // (2, 7), (4, 9) and (17, 29) are constrained
    assert_msg(num_bipartitions == expected_num_bipartitions, "Column", "Number of bipartition from column 2.");
}

void test_constrained_position_map_column2(Column* column) {
    const std::unordered_map<uint32_t, uint32_t>* map = column->get_constrained_position_map();

    /**
     * Clusters are {2, 3, 4, 7, 8, 9, 17, 29}
     *  Constrained Pair (2, 7) -> Position of 7 is 3 and Position of 2 is 0. So map is 0 -> 3.
     *  Constrained Pair (4, 9) -> Position of 9 is 5 and Position of 4 is 2. So map is 2 -> 5.
     *  Constrained Pair (17, 29) -> Position of 29 is 7 and Position of 17 is 6. So map is 6 -> 7.
     */

    assert_msg(map->size() == 3, "Column", "Size of column 2 constrained position map.");

    assert(map->at(0) == 3);
    assert(map->at(2) == 5);
    assert(map->at(6) == 7);
    assert_msg(true, "Column", "Content of column 2 constrained position map.");
}

void test_sorted_free_read_cluster_positions_column2(Column* column) {
    const std::vector<uint32_t>* positions = column->get_sorted_free_read_cluster_positions();
    
    /**
     * Clusters are {2, 3, 4, 7, 8, 9, 17, 29}
     * Free Read Clusters are {2, 3, 4, 8, 17}
     * Number of Reads associated with these clusters are {2, 2, 5, 2, 3}
     * 
     * If we sort (in ascending order) our free read clusters based on how many reads they have we get {2, 3, 8, 17, 4}
     * Their positions in the original read clusters vector are {0, 1, 4, 6, 2}
     */
    
    std::vector<uint32_t> expected = {0, 1, 4, 6, 2};
    assert_msg(positions->size() == expected.size(), "Column", "Size of sorted free read cluster positions of column 2.");
    for (uint32_t i = 0; i < positions->size(); i++) {
        assert(positions->at(i) == expected[i]);
    }
    assert_msg(true, "Column", "Content of sorted free read cluster positions of column 2.");

}

void test_cluster_id_to_read_index_map_column2(Column* column) {
    const std::unordered_map<uint32_t, std::vector<uint32_t>>* cluster_map = column->get_cluster_id_to_read_index_map();
    
    /**
     * Expected mapping:
     *  Cluster 8 -> read indices {0, 2}
     *  Cluster 7 -> read indices {1}
     *  Cluster 17 -> read indices {3, 12}
     *  Cluster 3 -> read indices {4, 7}
     *  Cluster 9 -> read indices {5, 6, 10}
     *  Cluster 2 -> read indices {8}
     *  Cluster 4 -> read indices {9, 11}
     *  Cluster 29 -> read indices {13}
     */
    assert(cluster_map->size() == 8);
    
    assert(cluster_map->at(8).size() == 2);
    assert(cluster_map->at(8)[0] == 0);
    assert(cluster_map->at(8)[1] == 2);

    assert(cluster_map->at(7).size() == 1);
    assert(cluster_map->at(7)[0] == 1);

    assert(cluster_map->at(17).size() == 2);
    assert(cluster_map->at(17)[0] == 3);
    assert(cluster_map->at(17)[1] == 12);
    
    assert(cluster_map->at(3).size() == 2);
    assert(cluster_map->at(3)[0] == 4);
    assert(cluster_map->at(3)[1] == 7);

    assert(cluster_map->at(9).size() == 3);
    assert(cluster_map->at(9)[0] == 5);
    assert(cluster_map->at(9)[1] == 6);
    assert(cluster_map->at(9)[2] == 10);

    assert(cluster_map->at(2).size() == 1);
    assert(cluster_map->at(2)[0] == 8);

    assert(cluster_map->at(4).size() == 2);
    assert(cluster_map->at(4)[0] == 9);
    assert(cluster_map->at(4)[1] == 11);

    assert(cluster_map->at(29).size() == 1);
    assert(cluster_map->at(29)[0] == 13);

    assert_msg(true, "Column", "Cluster ID to Read IDs map test passed for column 2.");
}

void test_column() {
    
    ReadSet* read_set = mock_readset_1();
    Column* column;
    
    const std::vector<uint32_t> curr_read_ids = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27};
    const std::vector<uint32_t> next_read_ids = {12, 14, 16, 17, 19, 20, 21, 22, 24, 25, 26, 27, 28, 29};

    column = new Column(curr_read_ids, next_read_ids, read_set);

    test_read_cluster_ids_column1(column);
    test_constrained_position_map_column1(column);
    test_cluster_id_to_read_index_map_column1(column);
    test_sorted_free_read_cluster_positions_column1(column);
    test_num_bipartitions_column1(column);
    test_bipartition_compatibility_column1(column);
    delete column;

    column = new Column(next_read_ids, {}, read_set);

    test_read_cluster_ids_column2(column);
    test_constrained_position_map_column2(column);
    test_cluster_id_to_read_index_map_column2(column);
    test_sorted_free_read_cluster_positions_column2(column);
    test_num_bipartitions_column2(column);
    delete column;

    // deleting objects
    delete read_set;
    
}