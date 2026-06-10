#include "test_phasingcolumniterator.h"

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
        std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next(true);
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
        std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next(true);
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
        std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next(true);
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


void test_phasingcolumniterator_unselected_reads() {

    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_2();

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
            std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next(true);
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
    variant_info_table[1].active_alleles[2] = false; // setting allele 2 for SV at 200 as inactive
    variant_info_table[1].phasable = true;
    variant_info_table[3].active_alleles[3] = false; // setting allele 3 for SV at 400 as inactive
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
            std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next(true);
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
    variant_info_table[3].active_alleles[2] = false; // setting allele 2 for SV at 400 as inactive
    variant_info_table[3].phasable = true;
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
            std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next(true);
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


void test_phasingcolumniterator_gapped_reads() {
    
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_2();
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

    variant_info_table[3].active_alleles[0] = false; variant_info_table[3].active_alleles[1] = false; variant_info_table[3].phasable = true;
    variant_info_table[1].active_alleles[2] = false; variant_info_table[1].phasable = true;
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
            std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next(true);
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
            std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next(true);
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

void test_jumping_columns() {
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_2();

    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read01", 60, 0); read_set->add(read1); read1->setSelected(false);
    Read* read2 = new Read("read02", 60, 0); read_set->add(read2); read2->setSelected(true);
    Read* read3 = new Read("read03", 60, 0); read_set->add(read3); read3->setSelected(true);
    Read* read4 = new Read("read04", 60, 0); read_set->add(read4); read4->setSelected(true);
    Read* read5 = new Read("read05", 60, 0); read_set->add(read5); read5->setSelected(true);
    Read* read6 = new Read("read06", 60, 0); read_set->add(read6); read6->setSelected(false);
    Read* read7 = new Read("read07", 60, 0); read_set->add(read7); read7->setSelected(true);
    Read* read8 = new Read("read08", 60, 0); read_set->add(read8); read8->setSelected(true);
    Read* read9 = new Read("read09", 60, 0); read_set->add(read9); read9->setSelected(true);
    Read* read10 = new Read("read10", 60, 0); read_set->add(read10); read10->setSelected(false);
    Read* read11 = new Read("read11", 60, 0); read_set->add(read11); read11->setSelected(true);
    Read* read12 = new Read("read12", 60, 0); read_set->add(read12); read12->setSelected(true);
    Read* read13 = new Read("read13", 60, 0); read_set->add(read13); read13->setSelected(true);

    /* Selected Reads */
    {
        read2->addVariant(100, std::vector<uint32_t>{10, 20});      // ALLELE1
        read2->addVariant(200, std::vector<uint32_t>{20, 5, 10});   // ALLELE1 -> allele at index 0 dropped
        
        read3->addVariant(100, std::vector<uint32_t>{10, 20});      // ALLELE1
        read3->addVariant(200, std::vector<uint32_t>{15, 8, 10});   // ALLELE1 -> allele at index 0 dropped

        read4->addVariant(100, std::vector<uint32_t>{10, 15});       // ALLELE1
        read4->addVariant(200, std::vector<uint32_t>{20, 5, 2});    // ALLELE2 -> allele at index 0 dropped

        read5->addVariant(100, std::vector<uint32_t>{10, 15});      // ALLELE1
        read5->addVariant(200, std::vector<uint32_t>{20, 5, 2});    // ALLELE2 -> allele at index 0 dropped

        read7->addVariant(200, std::vector<uint32_t>{15, 8, 2});    // ALLELE2 -> allele at index 0 dropped
        read7->addVariant(300, std::vector<uint32_t>{2, 5});        // ALLELE1
        
        read8->addVariant(200, std::vector<uint32_t>{15, 5, 8});    // ALLELE1 -> allele at index 0 dropped
        read8->addVariant(300, std::vector<uint32_t>{8, 5});        // ALLELE2
        read8->addVariant(400, std::vector<uint32_t>{5, 2, 90, 15});   // ALLELE2 -> allele at index 2, 3 dropped

        read9->addVariant(200, std::vector<uint32_t>{15, 8, 2});    // ALLELE2 -> allele at index 0 dropped
        read9->addVariant(300, std::vector<uint32_t>{2, 5});        // ALLELE1
        read9->addVariant(400, std::vector<uint32_t>{5, 2, 90, 15});    // ALLELE2 -> allele at index 2, 3 dropped
        
        read11->addVariant(400, std::vector<uint32_t>{5, 2, 90, 15});    // ALLELE2 -> allele at index 2, 3 dropped
        read11->addVariant(500, std::vector<uint32_t>{2, 5});        // ALLELE 1

        read12->addVariant(400, std::vector<uint32_t>{5, 2, 90, 15});    // ALLELE2 -> allele at index 2, 3 dropped
        read12->addVariant(500, std::vector<uint32_t>{5, 3});        // ALLELE 2

        read13->addVariant(400, std::vector<uint32_t>{5, 2, 90, 15});    // ALLELE2 -> allele at index 2, 3 dropped
        read13->addVariant(500, std::vector<uint32_t>{5, 3});        // ALLELE 2
        /**
         * Final superreads:
         * SR0 -> A1 - A1 - A2 - A2 - A2
         * SR1 -> A1 - A2 - A1 - A2 - A1
         */

    }
    /* Unselected Reads */
    {
        read1->addVariant(100, std::vector<uint32_t>{10, 20});      // ALLELE1

        read6->addVariant(200, std::vector<uint32_t>{15, 8, 10});   // ALLELE1

        read10->addVariant(300, std::vector<uint32_t>{15, 8});      // ALLELE2
        read10->addVariant(400, std::vector<uint32_t>{5, 2, 90, 15});   // ALLELE2 -> allele at index 2, 3 dropped
    }

    read_set->initialize();
    /* Checking if the IDs are correctly assigned */
    {
        assert(read1->getID() == 0);
        assert(read2->getID() == 1);
        assert(read3->getID() == 2);
        assert(read4->getID() == 3);
        assert(read5->getID() == 4);
        assert(read6->getID() == 5);
        assert(read7->getID() == 6);
        assert(read8->getID() == 7);
        assert(read9->getID() == 8);
        assert(read10->getID() == 9);
        assert(read11->getID() == 10);
        assert(read12->getID() == 11);
        assert(read13->getID() == 12);
    }
    /* Checking if the Position to Entry map is correct */
    {
        assert(read_set->TEST_get_pos_to_entry_map(100).size() ==  5);
        assert(read_set->TEST_get_pos_to_entry_map(200).size() ==  8);
        assert(read_set->TEST_get_pos_to_entry_map(300).size() ==  4);
        assert(read_set->TEST_get_pos_to_entry_map(400).size() ==  6);
        assert(read_set->TEST_get_pos_to_entry_map(500).size() ==  3);
    }
    for (size_t i = 0; i < variant_info_table.size(); ++i) {
        if (variant_info_table[i].phasable) {
			read_set->setEntryAlleles(variant_info_table[i].position, variant_info_table[i].active_alleles);
		}
    }
    /**
     * Setting up the conditions where the test for PhasingDPTable failed.
     * 
     * PhasingColumnIterator failed when jump_to_column(1) was executed after jump_to_column(2)
     * There should be no reads already in the active reads vector (so cannot execute after jump_to_column(0))
     *  so that the reads are all gathered by jump_to_column().
     */
    variant_info_table[1].genotype_likelihoods.increment_by_index(2, 0.5L); // 2 -> 1/1
    variant_info_table[1].genotype_likelihoods.increment_by_index(4, 0.3L); // 4 -> 1/2
    variant_info_table[1].genotype_likelihoods.increment_by_index(5, 0.2L); // 5 -> 2/2
    std::vector<uint32_t> selected_genotype_indices = variant_info_table[1].genotype_likelihoods.select_genotypes();
	variant_info_table[1].update_active_alleles(2, selected_genotype_indices);
    read_set->setEntryAlleles(variant_info_table[1].position, variant_info_table[1].active_alleles);
    read_set->resetTags();

    PhasingColumnIterator* iterator = new PhasingColumnIterator(*read_set, &variant_info_table, false);

    iterator->jump_to_column(1);
    std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next(true);
    assert(column->size() == 7);

    assert_msg(true, "PhasingColumnIterator", "Jumping columns.");
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
    variant_info_table[1].active_alleles[2] = false; // setting allele 2 for SV at 200 as inactive
    variant_info_table[1].phasable = true;
    variant_info_table[3].active_alleles[3] = false; // setting allele 3 for SV at 400 as inactive
    column_iterator = new PhasingColumnIterator(*read_set, &variant_info_table, false);
    second_round_phasing_tests(column_iterator);
    delete column_iterator;

    
    variant_info_table[3].active_alleles[2] = false; // setting allele 2 for SV at 400 as inactive
    variant_info_table[3].phasable = true;
    column_iterator = new PhasingColumnIterator(*read_set, &variant_info_table, false);
    third_round_phasing_tests(column_iterator);
    delete read_set;
    delete column_iterator;

    test_phasingcolumniterator_unselected_reads();

    test_phasingcolumniterator_gapped_reads();

    test_jumping_columns();
}