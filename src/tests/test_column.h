#ifndef TEST_COLUMN_H
#define TEST_COLUMN_H

/**
 * Unit tests for Column class.
 */

#include "../column.h"
#include "../tests_data.h"
#include <cassert>

void test_read_cluster_ids_column1(Column* column);

void test_num_bipartitions_column1(Column* column);

void test_constrained_position_map_column1(Column* column);

void test_sorted_free_read_cluster_positions_column1(Column* column);

void test_cluster_id_to_read_index_map_column1(Column* column);

void test_bipartition_compatibility_column1(Column* column);

void test_read_cluster_ids_column2(Column* column);

void test_num_bipartitions_column2(Column* column);

void test_constrained_position_map_column2(Column* column);

void test_sorted_free_read_cluster_positions_column2(Column* column);

void test_cluster_id_to_read_index_map_column2(Column* column);

void test_column();

#endif // TEST_COLUMN_H