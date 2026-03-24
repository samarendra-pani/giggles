#ifndef TEST_HAPLOTAGCOMPUTER_H
#define TEST_HAPLOTAGCOMPUTER_H

#include "../haplotagcomputer.h"
#include "../phasesetcomputer.h"
#include "../../../tests_data.h"
#include <cassert>

/** Position to index map testing. Just a sanity check. */
void test_position_to_index();

/** Reads that are tagged based on the partition from the DP table. Just a sanity check. */
void test_haplotag_selected_reads();

/** Calculation of edit distance from superreads. */
void test_distance_calculation_from_superreads();

/** Reads that were not used for phasing and hence don't have bipartition. */
void test_haplotag_unselected_reads();

void test_haplotagcomputer();

#endif // TEST_HAPLOTAGCOMPUTER_H