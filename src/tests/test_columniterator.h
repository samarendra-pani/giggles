#ifndef TEST_COLUMNITERATOR_H
#define TEST_COLUMNITERATOR_H

#include "../columniterator.h"
#include "../tests_data.h"
#include <cassert>


void test_mock_readset1() {

    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_1();
    ReadSet* read_set = mock_readset_1();

    ColumnIterator* column_iterator = new ColumnIterator(*read_set, &variant_info_table);
    
    assert_msg(column_iterator->get_column_count() == 3, "ColumnIterator", "Column count should be 3.");
    assert_msg(column_iterator->get_read_count() == read_set->size(), "ColumnIterator", "Read count should match read set size.");
    const std::vector<uint32_t>* positions = column_iterator->get_positions();
    assert_msg(positions->size() == 3, "ColumnIterator", "Positions size should be 3.");
    assert(positions->at(0) == 1);
    assert(positions->at(1) == 100);
    assert(positions->at(2) == 200);
    assert_msg(true, "ColumnIterator", "Correct positions.");

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
        assert_msg(true, "ColumnIterator", "Advanced to next column successfully.");
    }
    assert_msg(count == 3, "ColumnIterator", "Total columns iterated should be 3.");

    delete read_set;
    delete column_iterator;
}

void test_empty_columns() {
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_2();
    
    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3);
    std::vector<uint32_t> scores_1 = std::vector<uint32_t>{10, 90};
    read1->addVariant(100, scores_1);
    read2->addVariant(300, scores_1);
    read3->addVariant(500, scores_1);

    read_set->sort();
    read_set->reassignReadIds();

    assert(read1->getID() == 0);
    assert(read2->getID() == 1);
    assert(read3->getID() == 2);

    ColumnIterator* column_iterator = new ColumnIterator(*read_set, &variant_info_table);
    
    assert_msg(column_iterator->get_column_count() == 5, "ColumnIterator", "Column count for empty column example should be 5.");
    assert_msg(column_iterator->get_read_count() == read_set->size(), "ColumnIterator", "Read count for empty column example should match read set size.");
    const std::vector<uint32_t>* positions = column_iterator->get_positions();
    assert_msg(positions->size() == 5, "ColumnIterator", "Positions size for empty column example should be 5.");
    assert(positions->at(0) == 100);
    assert(positions->at(1) == 200);
    assert(positions->at(2) == 300);
    assert(positions->at(3) == 400);
    assert(positions->at(4) == 500);
    assert_msg(true, "ColumnIterator", "Correct positions for empty column example.");
    
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
        assert_msg(true, "ColumnIterator", "Advanced to next column for empty column example successfully.");
    }
    assert_msg(count == 5, "ColumnIterator", "Total columns iterated for empty column example should be 5.");

    delete read_set;
    delete column_iterator;
}

void test_gapped_reads() {
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_2();
    
    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3);
    std::vector<uint32_t> scores_1 = std::vector<uint32_t>{10, 90};
    read1->addVariant(100, scores_1); read1->addVariant(300, scores_1);
    read2->addVariant(300, scores_1);
    read3->addVariant(500, scores_1);

    read_set->sort();
    read_set->reassignReadIds();

    assert(read1->getID() == 0);
    assert(read2->getID() == 1);
    assert(read3->getID() == 2);

    ColumnIterator* column_iterator = new ColumnIterator(*read_set, &variant_info_table);
    
    assert_msg(column_iterator->get_column_count() == 5, "ColumnIterator", "Column count for gapped-read example should be 5.");
    assert_msg(column_iterator->get_read_count() == read_set->size(), "ColumnIterator", "Read count for gapped-read example should match read set size.");
    const std::vector<uint32_t>* positions = column_iterator->get_positions();
    assert_msg(positions->size() == 5, "ColumnIterator", "Positions size for gapped-read example should be 5.");
    assert(positions->at(0) == 100);
    assert(positions->at(1) == 200);
    assert(positions->at(2) == 300);
    assert(positions->at(3) == 400);
    assert(positions->at(4) == 500);
    assert_msg(true, "ColumnIterator", "Correct positions for gapped-read example.");
    
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
        assert_msg(true, "ColumnIterator", "Advanced to next column for gapped-read example successfully.");
    }
    assert_msg(count == 5, "ColumnIterator", "Total columns iterated for gapped-read example should be 5.");

    delete read_set;
    delete column_iterator;
}

void test_columniterator() {
    
    test_mock_readset1();
    test_empty_columns();
    test_gapped_reads();

}

#endif // TEST_COLUMNITERATOR_H