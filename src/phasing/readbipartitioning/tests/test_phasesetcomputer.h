#ifndef TEST_PHASESETCOMPUTER_H
#define TEST_PHASESETCOMPUTER_H

#include "../phasesetcomputer.h"
#include "../../../tests_data.h"
#include <cassert>

/** 
 * Testing the heterozygous positions determined.
 * NOTE: The test can not call a function. So we have copy-pasted part of compute_phasesets()
 */
void test_heterozygous_positions();

/** Reads do not overlap on a heterozygous position. Phaseset should remain same. */
void test_non_overlapping_reads();

/** Reads overlap on a heterozygous position. Phaseset should expand. */
void test_overlapping_reads();

/** Reads that do not cover any heterozygous position. */
void test_read_not_covering_het();

/** Reads that were not selected for phasing. */
void test_unselected_read();

void test_phasesetcomputer();

#endif // TEST_PHASESETCOMPUTER_H