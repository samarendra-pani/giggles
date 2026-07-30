#ifndef TEST_HAPLOTYPESAMPLER_H
#define TEST_HAPLOTYPESAMPLER_H

#include "../haplotypesampler.h"
#include "../../tests_data.h"
#include <cassert>

/* Testing haplotype sample for three columns. */
void test_multiplecolumns();

/* Testing haplotype sample for a single column (emission cost based). */
void test_singlecolumn();

void test_haplotypesampler();

#endif // TEST_HAPLOTYPESAMPLER_H