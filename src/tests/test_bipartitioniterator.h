#ifndef TEST_BIPARTITIONITERATOR_H
#define TEST_BIPARTITIONITERATOR_H

/**
 * Unit tests for the Biparition Iterator
 * This builds on top of the Gray Code tests.
 */

#include "../bipartitioniterator.h"
#include "../tests_data.h"
#include <cassert>

void test_constrained_bit_column1(BipartitionIterator* iterator);

void test_advancement_column1(BipartitionIterator* iterator);

void test_constrained_bit_column2(BipartitionIterator* iterator);

void test_advancement_column2(BipartitionIterator* iterator);

void test_bipartitioniterator();

#endif // TEST_BIPARTITIONITERATOR_H