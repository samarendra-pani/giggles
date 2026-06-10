#ifndef TEST_PHASINGDPTABLE_H
#define TEST_PHASINGDPTABLE_H

#include "../phasingdptable.h"
#include "../../tests_data.h"
#include <cassert>

/**
 * Single variant position. All reads at the position have been selected for phasing.
 */
void test_singleposition_allreadsselected();

/**
 * Single variant position. Some reads at the position have been selected for phasing.
 */
void test_singleposition_somereadsselected();

/**
 * Single SV variant position. Checking for when the position becomes considered for phasing. 
 */
void test_singleposition_sv();

/**
 * Multi-position variant read phasing with SVs.
 */
void test_multiposition();

void test_phasingdptable();

#endif // TEST_PHASINGDPTABLE_H