#ifndef TEST_ENTRY_H
#define TEST_ENTRY_H

/**
 * Unit tests for Entry class.
 */

#include "../entry.h"
#include "../tests_data.h"
#include <cassert>

void test_entry() {

    Entry* entry;
    entry = new Entry();
    assert_msg(entry->get_allele_type() == Entry::BLANK, "Entry", "Initialization allele type.");
    delete entry;

    entry = new Entry(1, std::vector<uint32_t>{10, 90});
    std::vector<long double> emission_scores = entry->get_emission_scores();
    assert_msg(emission_scores.size() == 2, "Entry", "Emission score size");
    assert_msg(emission_scores[0] < emission_scores[1], "Entry", "Relative emission scores.");
}

#endif //TEST_ENTRY_H