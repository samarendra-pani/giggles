#include "test_entry.h"

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

void test_entry() {

    test_entry_allele_type();
    test_entry_scores();
}