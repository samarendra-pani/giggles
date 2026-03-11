#ifndef TEST_PHASINGCOLUMNITERATOR_H
#define TEST_PHASINGCOLUMNITERATOR_H

#include "../phasingcolumniterator.h"
#include "../../tests_data.h"
#include <cassert>

void first_round_phasing_tests(PhasingColumnIterator* iterator);

void second_round_phasing_tests(PhasingColumnIterator* iterator);

void third_round_phasing_tests(PhasingColumnIterator* iterator);

void test_phasingcolumniterator_unselected_reads(std::vector<variant_information_t> variant_info_table);

void test_phasingcolumniterator_gapped_reads(std::vector<variant_information_t> variant_info_table);

void test_phasingcolumniterator();

#endif // TEST_PHASINGCOLUMNITERATOR_H