#ifndef TEST_VARIANTINFO_H
#define TEST_VARIANTINFO_H

/**
 * Unit tests for variant_information_t.
 */

#include "../tests_data.h"
#include <cassert>

void test_empty_variantinfo();

void test_active_allele_functions();

void test_genotype_likelihoods();

void test_update_active_allele();

void test_variantinfo();

#endif // TEST_VARIANTINFO_H