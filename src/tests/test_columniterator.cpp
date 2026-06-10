#include "test_columniterator.h"

void test_iterator_readset1() {

    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_1();
    ReadSet* read_set = mock_readset_1();

    ColumnIterator* column_iterator = new ColumnIterator(*read_set, &variant_info_table);
    
    assert_msg(column_iterator->get_column_count() == 3, "ColumnIterator", "Forward Iteration: Column count should be 3.");
    assert_msg(column_iterator->get_read_count() == read_set->size(), "ColumnIterator", "Forward Iteration: Read count should match read set size.");

    std::vector<std::vector<uint32_t>> expected_entry_id = {
        {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
        {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27},
        {12, 14, 16, 17, 19, 20, 21, 22, 24, 25, 26, 27, 28, 29}
    };
    
    uint32_t count = 0;
    while (column_iterator->has_next()) {
        std::unique_ptr<std::vector<const Entry*> > column = column_iterator->get_next();
        // check expected entries
        const std::vector<uint32_t>& expected_ids = expected_entry_id[count];
        assert(column->size() == expected_ids.size());
        for (size_t i = 0; i < expected_ids.size(); i++) {
            const Entry* entry = column->at(i);
            assert(entry->get_read_id() == expected_ids[i]);
        }
        count++;
        assert_msg(true, "ColumnIterator", "Forward Iteration: Advanced to next column successfully.");
    }
    assert_msg(count == 3, "ColumnIterator", "Forward Iteration: Total columns iterated should be 3.");

    column_iterator->jump_to_column(0);
    std::unique_ptr<std::vector<const Entry*> > column0 = column_iterator->get_next();
    const std::vector<uint32_t>& expected_ids0 = expected_entry_id[0];
    assert(column0->size() == expected_ids0.size());
    for (size_t i = 0; i < expected_ids0.size(); i++) {
        const Entry* entry = column0->at(i);
        assert(entry->get_read_id() == expected_ids0[i]);
    }
    assert_msg(true, "ColumnIterator", "Forward Iteration: Jump to Column 0.");

    column_iterator->jump_to_column(2);
    std::unique_ptr<std::vector<const Entry*> > column2 = column_iterator->get_next();
    const std::vector<uint32_t>& expected_ids2 = expected_entry_id[2];
    assert(column2->size() == expected_ids2.size());
    for (size_t i = 0; i < expected_ids2.size(); i++) {
        const Entry* entry = column2->at(i);
        assert(entry->get_read_id() == expected_ids2[i]);
    }
    assert_msg(true, "ColumnIterator", "Forward Iteration: Jump to Column 2.");

    column_iterator->jump_to_column(1);
    std::unique_ptr<std::vector<const Entry*> > column1 = column_iterator->get_next();
    const std::vector<uint32_t>& expected_ids1 = expected_entry_id[1];
    assert(column1->size() == expected_ids1.size());
    for (size_t i = 0; i < expected_ids1.size(); i++) {
        const Entry* entry = column1->at(i);
        assert(entry->get_read_id() == expected_ids1[i]);
    }
    assert_msg(true, "ColumnIterator", "Forward Iteration: Jump to Column 1.");

    delete read_set;
    delete column_iterator;
}

void test_iterator_empty_columns() {
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_2();
    
    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3);
    std::vector<uint32_t> scores_1 = std::vector<uint32_t>{10, 90};
    read1->addVariant(100, scores_1);
    read2->addVariant(300, scores_1);
    read3->addVariant(500, scores_1);

    read_set->initialize();
    
    assert(read1->getID() == 0);
    assert(read2->getID() == 1);
    assert(read3->getID() == 2);

    ColumnIterator* column_iterator = new ColumnIterator(*read_set, &variant_info_table);
    
    assert_msg(column_iterator->get_column_count() == 5, "ColumnIterator", "Forward Iteration: Column count for empty column example should be 5.");
    assert_msg(column_iterator->get_read_count() == read_set->size(), "ColumnIterator", "Forward Iteration: Read count for empty column example should match read set size.");
    
    std::vector<std::vector<uint32_t>> expected_entry_id = {{0}, {}, {1}, {}, {2}};
    
    uint32_t count = 0;
    while (column_iterator->has_next()) {
        std::unique_ptr<std::vector<const Entry*> > column = column_iterator->get_next();
        // check expected entries
        const std::vector<uint32_t>& expected_ids = expected_entry_id[count];
        assert(column->size() == expected_ids.size());
        for (size_t i = 0; i < expected_ids.size(); i++) {
            const Entry* entry = column->at(i);
            assert(entry->get_read_id() == expected_ids[i]);
        }
        count++;
        assert_msg(true, "ColumnIterator", "Forward Iteration: Advanced to next column for empty column example successfully.");
    }
    assert_msg(count == 5, "ColumnIterator", "Forward Iteration: Total columns iterated for empty column example should be 5.");

    delete read_set;
    delete column_iterator;
}

void test_iterator_gapped_reads() {
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_2();
    
    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3);
    std::vector<uint32_t> scores_1 = std::vector<uint32_t>{10, 90};
    read1->addVariant(100, scores_1); read1->addVariant(300, scores_1);
    read2->addVariant(300, scores_1);
    read3->addVariant(500, scores_1);

    read_set->initialize();
    
    assert(read1->getID() == 0);
    assert(read2->getID() == 1);
    assert(read3->getID() == 2);

    ColumnIterator* column_iterator = new ColumnIterator(*read_set, &variant_info_table);
    
    assert_msg(column_iterator->get_column_count() == 5, "ColumnIterator", "Forward Iteration: Column count for gapped-read example should be 5.");
    assert_msg(column_iterator->get_read_count() == read_set->size(), "ColumnIterator", "Forward Iteration: Read count for gapped-read example should match read set size.");
    
    std::vector<std::vector<uint32_t>> expected_entry_id = {{0}, {0}, {0, 1}, {}, {2}};
    
    uint32_t count = 0;
    while (column_iterator->has_next()) {
        std::unique_ptr<std::vector<const Entry*> > column = column_iterator->get_next();
        // check expected entries
        const std::vector<uint32_t>& expected_ids = expected_entry_id[count];
        assert(column->size() == expected_ids.size());
        for (size_t i = 0; i < expected_ids.size(); i++) {
            const Entry* entry = column->at(i);
            assert(entry->get_read_id() == expected_ids[i]);
        }
        count++;
        assert_msg(true, "ColumnIterator", "Forward Iteration: Advanced to next column for gapped-read example successfully.");
    }
    assert_msg(count == 5, "ColumnIterator", "Forward Iteration: Total columns iterated for gapped-read example should be 5.");

    delete read_set;
    delete column_iterator;
}

void test_backwarditerator_readset1() {
    
    ColumnIterator* column_iterator;
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_1();
    ReadSet* read_set = mock_readset_1();

    column_iterator = new ColumnIterator(*read_set, &variant_info_table);
    
    assert_msg(column_iterator->get_column_count() == 3, "ColumnIterator", "Backward Iteration: Column count should be 3.");
    assert_msg(column_iterator->get_read_count() == read_set->size(), "ColumnIterator", "Backward Iteration: Read count should match read set size.");

    std::vector<std::vector<uint32_t>> expected_entry_id = {
        {12, 14, 16, 17, 19, 20, 21, 22, 24, 25, 26, 27, 28, 29},
        {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27},
        {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
    };

    uint32_t count = 0;
    while (column_iterator->has_prev()) {
        std::unique_ptr<std::vector<const Entry*> > column = column_iterator->get_prev();
        // check expected entries
        const std::vector<uint32_t>& expected_ids = expected_entry_id[count];
        assert(column->size() == expected_ids.size());
        for (size_t i = 0; i < expected_ids.size(); i++) {
            const Entry* entry = column->at(i);
            assert(entry->get_read_id() == expected_ids[i]);
        }
        count++;
        assert_msg(true, "ColumnIterator", "Backward Iteration: Advanced to next column successfully.");
    }
    assert_msg(count == 3, "ColumnIterator", "Backward Iteration: Total columns iterated should be 3.");

    column_iterator->jump_to_column(0);
    std::unique_ptr<std::vector<const Entry*> > column0 = column_iterator->get_prev();
    const std::vector<uint32_t>& expected_ids0 = expected_entry_id[2];
    assert(column0->size() == expected_ids0.size());
    for (size_t i = 0; i < expected_ids0.size(); i++) {
        const Entry* entry = column0->at(i);
        assert(entry->get_read_id() == expected_ids0[i]);
    }
    assert_msg(true, "ColumnIterator", "Backward Iteration: Jump to Column 0.");

    column_iterator->jump_to_column(2);
    std::unique_ptr<std::vector<const Entry*> > column2 = column_iterator->get_prev();
    const std::vector<uint32_t>& expected_ids2 = expected_entry_id[0];
    assert(column2->size() == expected_ids2.size());
    for (size_t i = 0; i < expected_ids2.size(); i++) {
        const Entry* entry = column2->at(i);
        assert(entry->get_read_id() == expected_ids2[i]);
    }
    assert_msg(true, "ColumnIterator", "Backward Iteration: Jump to Column 2.");

    column_iterator->jump_to_column(1);
    std::unique_ptr<std::vector<const Entry*> > column1 = column_iterator->get_prev();
    const std::vector<uint32_t>& expected_ids1 = expected_entry_id[1];
    assert(column1->size() == expected_ids1.size());
    for (size_t i = 0; i < expected_ids1.size(); i++) {
        const Entry* entry = column1->at(i);
        assert(entry->get_read_id() == expected_ids1[i]);
    }
    assert_msg(true, "ColumnIterator", "Backward Iteration: Jump to Column 1.");

    delete read_set;
    delete column_iterator;
}

void test_backwarditerator_empty_columns() {
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_2();
    
    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3);
    std::vector<uint32_t> scores_1 = std::vector<uint32_t>{10, 90};
    read1->addVariant(100, scores_1);
    read2->addVariant(300, scores_1);
    read3->addVariant(500, scores_1);

    read_set->initialize();
    
    assert(read1->getID() == 0);
    assert(read2->getID() == 1);
    assert(read3->getID() == 2);

    ColumnIterator* column_iterator = new ColumnIterator(*read_set, &variant_info_table);
    
    assert_msg(column_iterator->get_column_count() == 5, "ColumnIterator", "Backward Iteration: Column count for empty column example should be 5.");
    assert_msg(column_iterator->get_read_count() == read_set->size(), "ColumnIterator", "Backward Iteration: Read count for empty column example should match read set size.");

    std::vector<std::vector<uint32_t>> expected_entry_id = {{2}, {}, {1}, {}, {0}};
    
    uint32_t count = 0;
    while (column_iterator->has_prev()) {
        std::unique_ptr<std::vector<const Entry*> > column = column_iterator->get_prev();
        // check expected entries
        const std::vector<uint32_t>& expected_ids = expected_entry_id[count];
        assert(column->size() == expected_ids.size());
        for (size_t i = 0; i < expected_ids.size(); i++) {
            const Entry* entry = column->at(i);
            assert(entry->get_read_id() == expected_ids[i]);
        }
        count++;
        assert_msg(true, "ColumnIterator", "Backward Iteration: Advanced to next column for empty column example successfully.");
    }
    assert_msg(count == 5, "ColumnIterator", "Backward Iteration: Total columns iterated for empty column example should be 5.");

    delete read_set;
    delete column_iterator;
}

void test_backwarditerator_gapped_reads() {
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_2();
    
    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3);
    std::vector<uint32_t> scores_1 = std::vector<uint32_t>{10, 90};
    read1->addVariant(100, scores_1); read1->addVariant(300, scores_1);
    read2->addVariant(300, scores_1);
    read3->addVariant(500, scores_1);

    read_set->initialize();
    
    assert(read1->getID() == 0);
    assert(read2->getID() == 1);
    assert(read3->getID() == 2);

    ColumnIterator* column_iterator = new ColumnIterator(*read_set, &variant_info_table);
    
    assert_msg(column_iterator->get_column_count() == 5, "ColumnIterator", "Backward Iteration: Column count for gapped-read example should be 5.");
    assert_msg(column_iterator->get_read_count() == read_set->size(), "ColumnIterator", "Backward Iteration: Read count for gapped-read example should match read set size.");
    
    std::vector<std::vector<uint32_t>> expected_entry_id = {{2}, {}, {0, 1}, {0}, {0}};
    
    uint32_t count = 0;
    while (column_iterator->has_prev()) {
        std::unique_ptr<std::vector<const Entry*> > column = column_iterator->get_prev();
        // check expected entries
        const std::vector<uint32_t>& expected_ids = expected_entry_id[count];
        assert(column->size() == expected_ids.size());
        for (size_t i = 0; i < expected_ids.size(); i++) {
            const Entry* entry = column->at(i);
            assert(entry->get_read_id() == expected_ids[i]);
        }
        count++;
        assert_msg(true, "ColumnIterator", "Backward Iteration: Advanced to next column for gapped-read example successfully.");
    }
    assert_msg(count == 5, "ColumnIterator", "Backward Iteration: Total columns iterated for gapped-read example should be 5.");

    delete read_set;
    delete column_iterator;
}

void test_mixediterations() {
    ColumnIterator* column_iterator;
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_1();
    ReadSet* read_set = mock_readset_1();

    column_iterator = new ColumnIterator(*read_set, &variant_info_table);

    std::vector<std::vector<uint32_t>> expected_entry_id = {
        {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
        {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27},
        {12, 14, 16, 17, 19, 20, 21, 22, 24, 25, 26, 27, 28, 29}
    };

    column_iterator->jump_to_column(2);
    assert(column_iterator->has_prev());
    assert(column_iterator->has_next());
    std::unique_ptr<std::vector<const Entry*> > column;
    column = column_iterator->get_next();
    assert(!column_iterator->has_next());
    assert(column->size() == expected_entry_id[2].size());
    for (size_t i = 0; i < expected_entry_id[2].size(); i++) {
        const Entry* entry = column->at(i);
        assert(entry->get_read_id() == expected_entry_id[2][i]);
    }

    column_iterator->jump_to_column(0);
    assert(column_iterator->has_prev());
    assert(column_iterator->has_next());
    column = column_iterator->get_prev();
    assert(!column_iterator->has_prev());
    assert(column->size() == expected_entry_id[0].size());
    for (size_t i = 0; i < expected_entry_id[0].size(); i++) {
        const Entry* entry = column->at(i);
        assert(entry->get_read_id() == expected_entry_id[0][i]);
    }
    
    column_iterator->jump_to_column(1);
    assert(column_iterator->has_prev());
    assert(column_iterator->has_next());
    column = column_iterator->get_prev();
    assert(column->size() == expected_entry_id[1].size());
    for (size_t i = 0; i < expected_entry_id[1].size(); i++) {
        const Entry* entry = column->at(i);
        assert(entry->get_read_id() == expected_entry_id[1][i]);
    }
    column = column_iterator->get_next();
    assert(column->size() == expected_entry_id[1].size());
    for (size_t i = 0; i < expected_entry_id[1].size(); i++) {
        const Entry* entry = column->at(i);
        assert(entry->get_read_id() == expected_entry_id[1][i]);
    }
    column = column_iterator->get_prev();
    assert(column->size() == expected_entry_id[0].size());
    for (size_t i = 0; i < expected_entry_id[0].size(); i++) {
        const Entry* entry = column->at(i);
        assert(entry->get_read_id() == expected_entry_id[0][i]);
    }
    column = column_iterator->get_next();
    assert(column->size() == expected_entry_id[2].size());
    for (size_t i = 0; i < expected_entry_id[2].size(); i++) {
        const Entry* entry = column->at(i);
        assert(entry->get_read_id() == expected_entry_id[2][i]);
    }
    assert(!column_iterator->has_prev());
    assert(!column_iterator->has_next());
    
    assert_msg(true, "ColumnIterator", "Column jumping and independence of get_next() and get_prev().");

}

void test_columniterator() {
    
    test_iterator_readset1();
    test_iterator_empty_columns();
    test_iterator_gapped_reads();

    test_backwarditerator_readset1();
    test_backwarditerator_empty_columns();
    test_backwarditerator_gapped_reads();

    test_mixediterations();

}