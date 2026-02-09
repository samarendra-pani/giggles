/**
 * TODO: New Assertion Macros
 */

#ifndef TEST_PHASINGCOLUMNITERATOR_H
#define TEST_PHASINGCOLUMNITERATOR_H

#include "../phasingcolumniterator.h"
#include "../../tests_data.h"
#include <cassert>


void first_round_phasing_tests(PhasingColumnIterator* iterator) {
    
    std::cout << "[PhasingColumnIterator] Testing for first round of phasing..." << std::endl;
    std::cout << "[PhasingColumnIterator] Testing advancement..." << std::endl;


    std::vector<std::vector<uint32_t>> expected_entry_id = {
        {0, 1},
        {0, 1, 2},
        {0, 1, 2, 3},
        {0, 1, 2, 3},
        {1, 3}
    };
    std::vector<std::vector<Entry::allele_t>> expected_alleles = {
        {Entry::ALLELE1, Entry::ALLELE1},
        {Entry::BLANK, Entry::BLANK, Entry::BLANK},
        {Entry::ALLELE2, Entry::ALLELE2, Entry::ALLELE2, Entry::ALLELE2},
        {Entry::BLANK, Entry::BLANK, Entry::BLANK, Entry::BLANK},
        {Entry::ALLELE1, Entry::ALLELE1}
    };
    uint32_t count = 0;
    while (iterator->has_next()) {
        std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next();
        std::cout << "[PhasingColumnIterator] Testing column " << count << "..." << std::endl;
        // check expected entries
        const std::vector<uint32_t>& expected_ids = expected_entry_id[count];
        const std::vector<Entry::allele_t>& expected_allele_types = expected_alleles[count];
        assert(column->size() == expected_ids.size());
        for (size_t i = 0; i < expected_ids.size(); i++) {
            const Entry* entry = column->at(i);
            assert(entry->get_read_id() == expected_ids[i]);
            assert(entry->get_allele_type() == expected_allele_types[i]);
        }
        count++;
    }    
    assert(count == 5);
}


void second_round_phasing_tests(PhasingColumnIterator* iterator) {
    
    std::cout << "[PhasingColumnIterator] Testing for non-first round of phasing (with BLANK)..." << std::endl;
    std::cout << "[PhasingColumnIterator] Testing advancement..." << std::endl;


    std::vector<std::vector<uint32_t>> expected_entry_id = {
        {0, 1},
        {0, 1, 2},
        {0, 1, 2, 3},
        {0, 1, 2, 3},
        {1, 3}
    };
    std::vector<std::vector<Entry::allele_t>> expected_alleles = {
        {Entry::ALLELE1, Entry::ALLELE1},
        {Entry::ALLELE1, Entry::ALLELE1, Entry::ALLELE1},
        {Entry::ALLELE2, Entry::ALLELE2, Entry::ALLELE2, Entry::ALLELE2},
        {Entry::BLANK, Entry::BLANK, Entry::BLANK, Entry::BLANK},
        {Entry::ALLELE1, Entry::ALLELE1}
    };
    uint32_t count = 0;
    while (iterator->has_next()) {
        std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next();
        std::cout << "[PhasingColumnIterator] Testing column " << count << "..." << std::endl;
        // check expected entries
        const std::vector<uint32_t>& expected_ids = expected_entry_id[count];
        const std::vector<Entry::allele_t>& expected_allele_types = expected_alleles[count];
        assert(column->size() == expected_ids.size());
        for (size_t i = 0; i < expected_ids.size(); i++) {
            const Entry* entry = column->at(i);
            assert(entry->get_read_id() == expected_ids[i]);
            assert(entry->get_allele_type() == expected_allele_types[i]);
        }
        count++;
    }    
    assert(count == 5);
}

void third_round_phasing_tests(PhasingColumnIterator* iterator) {
    
    std::cout << "[PhasingColumnIterator] Testing for non-first round of phasing (with EQUAL SCORES)..." << std::endl;
    std::cout << "[PhasingColumnIterator] Testing advancement..." << std::endl;


    std::vector<std::vector<uint32_t>> expected_entry_id = {
        {0, 1},
        {0, 1, 2},
        {0, 1, 2, 3},
        {0, 1, 2, 3},
        {1, 3}
    };
    std::vector<std::vector<Entry::allele_t>> expected_alleles = {
        {Entry::ALLELE1, Entry::ALLELE1},
        {Entry::ALLELE1, Entry::ALLELE1, Entry::ALLELE1},
        {Entry::ALLELE2, Entry::ALLELE2, Entry::ALLELE2, Entry::ALLELE2},
        {Entry::EQUAL_SCORES, Entry::EQUAL_SCORES, Entry::EQUAL_SCORES, Entry::EQUAL_SCORES},
        {Entry::ALLELE1, Entry::ALLELE1}
    };
    uint32_t count = 0;
    while (iterator->has_next()) {
        std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next();
        std::cout << "[PhasingColumnIterator] Testing column " << count << "..." << std::endl;
        // check expected entries
        const std::vector<uint32_t>& expected_ids = expected_entry_id[count];
        const std::vector<Entry::allele_t>& expected_allele_types = expected_alleles[count];
        assert(column->size() == expected_ids.size());
        for (size_t i = 0; i < expected_ids.size(); i++) {
            const Entry* entry = column->at(i);
            assert(entry->get_read_id() == expected_ids[i]);
            assert(entry->get_allele_type() == expected_allele_types[i]);
        }
        count++;
    }    
    
}



void test_phasingcolumniterator() {
    
    
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_2();
    ReadSet* read_set = mock_readset_2();
    
    PhasingColumnIterator* column_iterator = new PhasingColumnIterator(*read_set, &variant_info_table, true);
    
    std::cout << "[PhasingColumnIterator] Testing get_column_count..." << std::endl;
    assert(column_iterator->get_column_count() == 2);
    std::cout << "[PhasingColumnIterator] Testing get_read_count..." << std::endl;
    assert(column_iterator->get_read_count() == read_set->size());
    std::cout << "[PhasingColumnIterator] Testing get_positions..." << std::endl;
    const std::vector<uint32_t>* positions = column_iterator->get_positions();
    assert(positions->at(0) == 100);
    assert(positions->at(1) == 200);
    assert(positions->at(0) == 300);
    assert(positions->at(1) == 400);
    assert(positions->at(0) == 500);
    
    first_round_phasing_tests(column_iterator);
    delete column_iterator;

    // setting the active alleles for SVs
    variant_info_table[1].set_allele_inactive(2); // setting variant 0 for SV at 200 as inactive
    variant_info_table[3].set_allele_inactive(3); // setting variant 3 for SV at 400 as inactive
    PhasingColumnIterator* column_iterator = new PhasingColumnIterator(*read_set, &variant_info_table, false);
    second_round_phasing_tests(column_iterator);
    delete column_iterator;

    
    variant_info_table[3].set_allele_inactive(2); // setting variant 1 for SV at 400 as inactive
    PhasingColumnIterator* column_iterator = new PhasingColumnIterator(*read_set, &variant_info_table, false);
    third_round_phasing_tests(column_iterator);

    
    delete read_set;
    delete column_iterator;

}

#endif // TEST_PHASINGCOLUMNITERATOR_H