#ifndef TEST_PHASINGCOLUMNCOSTCOMPUTER_H
#define TEST_PHASINGCOLUMNCOSTCOMPUTER_H

#include "../phasingcolumncostcomputer.h"
#include "../phasingcolumniterator.h"
#include "../../tests_data.h"
#include <cassert>


void test_setting_partition();

void test_update_partition();

void test_equal_scores();

void test_multiallelic_site();

/**
 * TODO: More tests for cases where only haplotype has EQUAL SCORES
 */
void test_phasingcolumncostcomputer();

#endif // TEST_PHASINGCOLUMNCOSTCOMPUTER_H
