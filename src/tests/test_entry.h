#ifndef TEST_ENTRY_H
#define TEST_ENTRY_H

/**
 * Unit tests for Entry class.
 */

#include "../entry.h"
#include "../tests_data.h"
#include <cassert>

void test_entry_allele_type() {
    
    Entry* entry;
    entry = new Entry();
    assert_msg(entry->get_allele_type() == Entry::BLANK, "Entry", "Empty nitialization allele type.");
    delete entry;

    entry = new Entry(1, std::vector<uint32_t>{10, 90});
    assert_msg(entry->get_allele_type() == Entry::BLANK, "Entry", "Non-empty initialization allele type.");
    entry->set_allele_type({true, true});
    assert_msg(entry->get_allele_type() == Entry::ALLELE1, "Entry", "Setting allele type 1.");
    delete entry;
    
    entry = new Entry(1, std::vector<uint32_t>{90, 10});
    entry->set_allele_type({true, true});
    assert_msg(entry->get_allele_type() == Entry::ALLELE2, "Entry", "Setting allele type 2.");
    delete entry;

    entry = new Entry(1, std::vector<uint32_t>{50, 500});
    assert(entry->get_allele_type() == Entry::BLANK);
    entry->set_allele_type({true, true});
    assert(entry->get_allele_type() == Entry::EQUAL_SCORES);
    delete entry;
    entry = new Entry(1, std::vector<uint32_t>{50, 50});
    assert(entry->get_allele_type() == Entry::BLANK);
    entry->set_allele_type({true, true});
    assert(entry->get_allele_type() == Entry::EQUAL_SCORES);
    delete entry;
    entry = new Entry(1, std::vector<uint32_t>{500, 70});
    assert(entry->get_allele_type() == Entry::BLANK);
    entry->set_allele_type({true, true});
    assert(entry->get_allele_type() == Entry::EQUAL_SCORES);
    delete entry;
    assert_msg(true, "Entry", "Setting allele type EQUAL_SCORES.");
    

    entry = new Entry(1, std::vector<uint32_t>{10, 90, 5, 20});
    assert(entry->get_allele_type() == Entry::BLANK);
    entry->set_allele_type({true, false, true, false});
    assert_msg(entry->get_allele_type() == Entry::ALLELE2, "Entry", "Setting allele type with more than 2 alleles");
    delete entry;

    
}

void test_entry_scores() {

    Entry* entry;
    entry = new Entry(1, std::vector<uint32_t>{10, 90});
    std::vector<long double> emission_scores = entry->get_emission_scores();
    assert_msg(emission_scores.size() == 2, "Entry", "Emission score size");
    assert_msg(emission_scores[0] > emission_scores[1], "Entry", "Relative emission scores.");
    delete entry;

}

void test_entry_genotypes() {
    Entry* entry;
    
    entry = new Entry(1, std::vector<uint32_t>{10, 90});
    entry->set_allele_type({true, true});
    assert(entry->get_allele_type() == Entry::ALLELE1);
    assert_msg(entry->get_allele() == 0, "Entry", "Getting allele 0 as ALLELE1.");
    delete entry;

    entry = new Entry(1, std::vector<uint32_t>{90, 10});
    entry->set_allele_type({true, true});
    assert(entry->get_allele_type() == Entry::ALLELE2);
    assert_msg(entry->get_allele() == 1, "Entry", "Getting allele 1 as ALLELE2.");
    delete entry;

    entry = new Entry(1, std::vector<uint32_t>{50, 50});
    entry->set_allele_type({true, true});
    assert(entry->get_allele_type() == Entry::EQUAL_SCORES);
    assert_msg(entry->get_allele() == (uint32_t)-1, "Entry", "Getting (uint32_t)-1 for EQUAL_SCORES.");
    delete entry;

    entry = new Entry(1, std::vector<uint32_t>{10, 90, 5, 7});
    entry->set_allele_type({false, false, true, true});
    assert(entry->get_allele_type() == Entry::ALLELE1);
    assert(entry->get_allele() == 2);
    entry->set_allele_type({true, false, false, true});
    assert(entry->get_allele_type() == Entry::ALLELE2);
    assert(entry->get_allele() == 3);
    delete entry;
    entry = new Entry(1, std::vector<uint32_t>{50, 90, 5, 7});
    entry->set_allele_type({true, false, false, true});
    assert(entry->get_allele_type() == Entry::ALLELE2);
    assert(entry->get_allele() == 3);
    entry->set_allele_type({true, true, false, false});
    assert(entry->get_allele_type() == Entry::EQUAL_SCORES);
    assert(entry->get_allele() == (uint32_t)-1);
    delete entry;
    assert_msg(true, "Entry", "Alleles from genotypes in Entries with more than 2 alleles.");
    
}

void test_entry() {

    test_entry_allele_type();
    test_entry_scores();
    test_entry_genotypes();
    
}

#endif //TEST_ENTRY_H