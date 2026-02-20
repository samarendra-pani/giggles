#ifndef TEST_ENTRY_H
#define TEST_ENTRY_H

/**
 * Unit tests for Entry class.
 */

#include "../entry.h"
#include <cassert>

void test_entry() {

    Entry* entry = new Entry();
    assert(entry->get_allele_type() == Entry::BLANK);
    delete entry;

    Entry* entry = new Entry(1, std::vector<uint32_t>{10, 90});
    std::vector<long double> emission_scores = entry->get_emission_scores();
    assert(emission_scores.size() == 2);
    assert(emission_scores[0] < emission_scores[1]);
}

#endif //TEST_ENTRY_H