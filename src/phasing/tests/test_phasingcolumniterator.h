#ifndef TEST_PHASINGCOLUMNITERATOR_H
#define TEST_PHASINGCOLUMNITERATOR_H

#include "../phasingcolumniterator.h"
#include "../../tests_data.h"
#include <cassert>


void first_round_phasing_tests(PhasingColumnIterator* iterator) {
    
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
        {Entry::EQUAL_SCORES, Entry::EQUAL_SCORES}
    };
    uint32_t count = 0;
    while (iterator->has_next()) {
        std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next();
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
    assert_msg(true, "PhasingColumnIterator", "Iterator advancement for first round of phasing.");
    assert_msg(count == 5, "PhasingColumnIterator", "Advanced 5 times.");
}


void second_round_phasing_tests(PhasingColumnIterator* iterator) {
    
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
        {Entry::EQUAL_SCORES, Entry::EQUAL_SCORES}
    };
    uint32_t count = 0;
    while (iterator->has_next()) {
        std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next();
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
    assert_msg(true, "PhasingColumnIterator", "Iterator advancement for second round of phasing.");
    assert_msg(count == 5, "PhasingColumnIterator", "Advanced 5 times.");
}

void third_round_phasing_tests(PhasingColumnIterator* iterator) {
    
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
        {Entry::ALLELE1, Entry::ALLELE1, Entry::ALLELE1, Entry::ALLELE1},
        {Entry::EQUAL_SCORES, Entry::EQUAL_SCORES}
    };
    uint32_t count = 0;
    while (iterator->has_next()) {
        std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next();
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
    assert_msg(true, "PhasingColumnIterator", "Iterator advancement for third round of phasing.");
    assert_msg(count == 5, "PhasingColumnIterator", "Advanced 5 times."); 
}


void test_phasingcolumniterator_unselected_reads(std::vector<variant_information_t> variant_info_table) {

    ReadSet* read_set = mock_readset_2();
    read_set->getByName("read2", 0)->setSelected(false);
    read_set->getByName("read3", 0)->setSelected(false);
    PhasingColumnIterator* iterator;
    iterator = new PhasingColumnIterator(*read_set, &variant_info_table, true);

    // First round of phasing
    {
        std::vector<std::vector<uint32_t>> expected_entry_id = {
            {0},
            {0},
            {0, 3},
            {0, 3},
            {3}
        };
        std::vector<std::vector<Entry::allele_t>> expected_alleles = {
            {Entry::ALLELE1},
            {Entry::BLANK},
            {Entry::ALLELE2, Entry::ALLELE2},
            {Entry::BLANK, Entry::BLANK},
            {Entry::EQUAL_SCORES}
        };
        uint32_t count = 0;
        while (iterator->has_next()) {
            std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next();
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
    delete iterator;

    // Second round of phasing
    read_set->getByName("read3", 0)->setSelected(true);
    variant_info_table[1].set_allele_inactive(2); // setting allele 2 for SV at 200 as inactive
    variant_info_table[3].set_allele_inactive(3); // setting allele 3 for SV at 400 as inactive
    iterator = new PhasingColumnIterator(*read_set, &variant_info_table, false);
    {
        std::vector<std::vector<uint32_t>> expected_entry_id = {
            {0},
            {0, 2},
            {0, 2, 3},
            {0, 2, 3},
            {3}
        };
        std::vector<std::vector<Entry::allele_t>> expected_alleles = {
            {Entry::ALLELE1},
            {Entry::ALLELE1, Entry::ALLELE1},
            {Entry::ALLELE2, Entry::ALLELE2, Entry::ALLELE2},
            {Entry::BLANK, Entry::BLANK, Entry::BLANK},
            {Entry::EQUAL_SCORES}
        };
        uint32_t count = 0;
        while (iterator->has_next()) {
            std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next();
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
    delete iterator;

    // third round of phasing
    read_set->getByName("read2", 0)->setSelected(true);
    variant_info_table[3].set_allele_inactive(2); // setting allele 2 for SV at 400 as inactive
    iterator = new PhasingColumnIterator(*read_set, &variant_info_table, false);
    {
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
            {Entry::ALLELE1, Entry::ALLELE1, Entry::ALLELE1, Entry::ALLELE1},
            {Entry::EQUAL_SCORES, Entry::EQUAL_SCORES}
        };
        uint32_t count = 0;
        while (iterator->has_next()) {
            std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next();
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
    assert_msg(true, "PhasingColumnIterator", "Iterator test with unselected reads.");
}


void test_phasingcolumniterator_gapped_reads(std::vector<variant_information_t> variant_info_table) {
    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(true);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(true);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(true);
    Read* read4 = new Read("read4", 60, 0); read_set->add(read4); read4->setSelected(true);

    /**
     * Adding variants to the reads
     */
    std::vector<uint32_t> scores_1 = std::vector<uint32_t>{10, 90};
    std::vector<uint32_t> scores_2 = std::vector<uint32_t>{20, 30, 50};
    std::vector<uint32_t> scores_3 = std::vector<uint32_t>{85, 15};
    std::vector<uint32_t> scores_4 = std::vector<uint32_t>{5, 25, 35, 35};
    std::vector<uint32_t> scores_5 = std::vector<uint32_t>{40, 60};
    
    read1->addVariant(100, scores_1); read1->addVariant(200, scores_2); read1->addVariant(400, scores_4);
    read2->addVariant(100, scores_1); read2->addVariant(200, scores_2); read2->addVariant(500, scores_5);
    read3->addVariant(200, scores_2); read3->addVariant(400, scores_4);
    read4->addVariant(300, scores_3); read4->addVariant(400, scores_4); read4->addVariant(500, scores_5);

    /**
     * Need to set this manually since hash function to break ties does the tie breaking in weird way
     */
    read1->setID(0);
    read2->setID(1);
    read3->setID(2);
    read4->setID(3);

    /**
     * Reads Summary:
     * ID | Name   | Variants                 
     * ---|--------|--------------------------
     * 0  | read1  | 100, 200, 400       
     * 1  | read2  | 100, 200, 500
     * 2  | read3  | 200, 400
     * 3  | read4  | 300, 400, 500
     */

    variant_info_table[3].set_allele_inactive(0); variant_info_table[3].set_allele_inactive(1);
    variant_info_table[1].set_allele_inactive(2);
    PhasingColumnIterator* iterator;
    iterator = new PhasingColumnIterator(*read_set, &variant_info_table, true);

    // first round of phasing
    {
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
            {Entry::BLANK, Entry::BLANK, Entry::BLANK, Entry::ALLELE2},
            {Entry::BLANK, Entry::BLANK, Entry::BLANK, Entry::BLANK},
            {Entry::EQUAL_SCORES, Entry::EQUAL_SCORES}
        };
        uint32_t count = 0;
        while (iterator->has_next()) {
            std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next();
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

    // second round of phasing
    {
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
            {Entry::BLANK, Entry::BLANK, Entry::BLANK, Entry::ALLELE2},
            {Entry::EQUAL_SCORES, Entry::BLANK, Entry::EQUAL_SCORES, Entry::EQUAL_SCORES},
            {Entry::EQUAL_SCORES, Entry::EQUAL_SCORES}
        };
        uint32_t count = 0;
        while (iterator->has_next()) {
            std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next();
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
    assert_msg(true, "PhasingColumnIterator", "Iterator test with gapped reads.");
}


void test_phasingcolumniterator() {
    
    
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_2();
    ReadSet* read_set = mock_readset_2();
    PhasingColumnIterator* column_iterator;

    column_iterator = new PhasingColumnIterator(*read_set, &variant_info_table, true);
    
    assert_msg(column_iterator->get_column_count() == 5, "PhasingColumnIterator", "Column count should be 5.");
    assert_msg(column_iterator->get_read_count() == read_set->size(), "PhasingColumnIterator", "Read count should match read set size.");
    assert(column_iterator->get_position(0) == 100);
    assert(column_iterator->get_position(1) == 200);
    assert(column_iterator->get_position(2) == 300);
    assert(column_iterator->get_position(3) == 400);
    assert(column_iterator->get_position(4) == 500);
    assert_msg(true, "PhasingColumnIterator", "Positions of the iterator.");

    first_round_phasing_tests(column_iterator);
    delete column_iterator;

    // setting the active alleles for SVs
    variant_info_table[1].set_allele_inactive(2); // setting allele 2 for SV at 200 as inactive
    variant_info_table[3].set_allele_inactive(3); // setting allele 3 for SV at 400 as inactive
    column_iterator = new PhasingColumnIterator(*read_set, &variant_info_table, false);
    second_round_phasing_tests(column_iterator);
    delete column_iterator;

    
    variant_info_table[3].set_allele_inactive(2); // setting allele 2 for SV at 400 as inactive
    column_iterator = new PhasingColumnIterator(*read_set, &variant_info_table, false);
    third_round_phasing_tests(column_iterator);
    delete read_set;
    delete column_iterator;

    variant_info_table = mock_variant_info_table_2();
    test_phasingcolumniterator_unselected_reads(variant_info_table);

    variant_info_table = mock_variant_info_table_2();
    test_phasingcolumniterator_gapped_reads(variant_info_table);
}

#endif // TEST_PHASINGCOLUMNITERATOR_H