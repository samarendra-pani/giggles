#include "test_entry.h"

void test_entry_allele_type() {
    
    Entry* entry;
    entry = new Entry();
    assert_msg(entry->get_allele_type() == Entry::BLANK, "Entry", "Empty initialization allele type.");
    delete entry;

    entry = new Entry(1, std::vector<float>{0.9, 0.1});
    assert_msg(entry->get_allele_type() == Entry::BLANK, "Entry", "Non-empty initialization allele type.");
    entry->set_allele_type({true, true});
    assert_msg(entry->get_allele_type() == Entry::ALLELE1, "Entry", "Setting allele type 1.");
    delete entry;
    
    entry = new Entry(1, std::vector<float>{0.1, 0.9});
    entry->set_allele_type({true, true});
    assert_msg(entry->get_allele_type() == Entry::ALLELE2, "Entry", "Setting allele type 2.");
    delete entry;

    entry = new Entry(1, std::vector<float>{0.0f, 0.0f});
    assert(entry->get_allele_type() == Entry::BLANK);
    entry->set_allele_type({true, true});
    assert(entry->get_allele_type() == Entry::EQUAL_SCORES);
    delete entry;
    entry = new Entry(1, std::vector<float>{0.9f, 0.9f});
    assert(entry->get_allele_type() == Entry::BLANK);
    entry->set_allele_type({true, true});
    assert(entry->get_allele_type() == Entry::EQUAL_SCORES);
    delete entry;
    assert_msg(true, "Entry", "Setting allele type EQUAL_SCORES.");
    

    entry = new Entry(1, std::vector<float>{0.9, 0.1, 0.95, 0.8});
    assert(entry->get_allele_type() == Entry::BLANK);
    entry->set_allele_type({true, false, true, false});
    assert_msg(entry->get_allele_type() == Entry::ALLELE2, "Entry", "Setting allele type with more than 2 alleles");
    delete entry;   
}

void test_entry_scores() {

    Entry* entry;
    Entry::initialize_probability_cache(30.0f);
    entry = new Entry(1, std::vector<float>{0.9, 0.1});
    assert_msg(entry->get_emission_score(0) > entry->get_emission_score(1), "Entry", "Relative emission scores.");
    delete entry;
}

void test_entry() {

    test_entry_allele_type();
    test_entry_scores();
}