#ifndef TEST_PHASINGCOLUMNITERATOR_H
#define TEST_PHASINGCOLUMNITERATOR_H

#include "../phasingcolumniterator.h"
#include "../../tests_data.h"
#include <cassert>

void first_round_phasing_tests(PhasingColumnIterator* iterator);

void second_round_phasing_tests(PhasingColumnIterator* iterator);

void third_round_phasing_tests(PhasingColumnIterator* iterator);

void test_phasingcolumniterator_unselected_reads();

void test_phasingcolumniterator_gapped_reads();

void test_jumping_columns();

void test_phasingcolumniterator();

#endif // TEST_PHASINGCOLUMNITERATOR_H