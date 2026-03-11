#ifndef TEST_GENOTYPE_H
#define TEST_GENOTYPE_H

#include "../genotype.h"
#include "../tests_data.h"

void test_empty_genotype();

void test_diploid_construction_and_ordering();

void test_canonical_index_mapping();

void test_haploid_chrY();

void test_predicates();

void test_limits_and_exceptions();

void test_comparison_operators();

void test_genotype();


#endif // TEST_GENOTYPE_H