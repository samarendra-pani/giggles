#ifndef TEST_GRAYCODES_H
#define TEST_GRAYCODES_H

/**
 * Unit tests for GrayCodes utility functions.
 */

#include "../graycodes.h"
#include "../tests_data.h"
#include <cassert>

void test_zero_length_graycode();

void test_l_length_graycode(int l);

void test_length_4_graycode();

void test_graycodes();


#endif // TEST_GRAYCODES_H