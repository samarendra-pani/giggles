#ifndef TEST_COLUMNITERATOR_H
#define TEST_COLUMNITERATOR_H

#include "../columniterator.h"
#include "../tests_data.h"
#include <cassert>

void test_iterator_readset1();

void test_iterator_empty_columns();

void test_iterator_gapped_reads();

void test_backwarditerator_readset1();

void test_backwarditerator_empty_columns();

void test_backwarditerator_gapped_reads();

void test_mixediterations();

void test_columniterator();

#endif // TEST_COLUMNITERATOR_H