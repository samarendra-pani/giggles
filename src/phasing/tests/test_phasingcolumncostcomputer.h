#ifndef TEST_PHASINGCOLUMNCOSTCOMPUTER_H
#define TEST_PHASINGCOLUMNCOSTCOMPUTER_H

#include "../phasingcolumncostcomputer.h"
#include "../phasingcolumniterator.h"
#include "../../tests_data.h"
#include <cassert>

/**
 * Unphasable variant position should give
 * EQUAL SCORES and 0 cost.
 */
void test_unphasable_position();

/**
 * Position determined to be homozygous by genotypelikelihoods
 * ALLELE1 and 0 cost.
 */
void test_homozygous_position();

/**
 * Setting of partitions with no genotype likelihoods.
 * This allows all genotypes possible.
 */
void test_setting_partition_no_gl();

/**
 * Setting of partitions with no genotype likelihoods.
 * This imposes genotype restrictions.
 */
void test_setting_partition_with_gl();

void test_update_partition_no_gl();

/**
 * Testing some extra cases of EQUAL_SCORES
 */
void test_equal_scores();

void test_phasingcolumncostcomputer();

#endif // TEST_PHASINGCOLUMNCOSTCOMPUTER_H
