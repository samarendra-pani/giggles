#ifndef TEST_COLUMNITERATOR_H
#define TEST_COLUMNITERATOR_H

#include "../columniterator.h"
#include "../tests_data.h"
#include <cassert>


void test_columniterator() {
    
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_1();
    ReadSet* read_set = mock_readset_1();

    ColumnIterator* column_iterator = new ColumnIterator(*read_set, &variant_info_table);
    
    assert_msg(column_iterator->get_column_count() == 2, "ColumnIterator", "Column count test.");
    assert_msg(column_iterator->get_read_count() == read_set->size(), "ColumnIterator", "Read count test.");
    const std::vector<uint32_t>* positions = column_iterator->get_positions();
    assert_msg(positions->size() == 2, "ColumnIterator", "Positions size test.");
    assert(positions->at(0) == 100);
    assert(positions->at(1) == 200);
    assert_msg(true, "ColumnIterator", "Correct positions");

    std::cout << "[ColumnIterator] Testing advancement..." << std::endl;
    
    std::vector<std::vector<uint32_t>> expected_entry_id = {
        {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27},
        {12, 14, 16, 17, 19, 20, 21, 22, 24, 25, 26, 27, 28, 29}
    };
    
    uint32_t count = 0;
    while (column_iterator->has_next()) {
        std::unique_ptr<std::vector<const Entry*> > column = column_iterator->get_next();
        std::cout << "[ColumnIterator] Testing column " << count << "..." << std::endl;
        // check expected entries
        const std::vector<uint32_t>& expected_ids = expected_entry_id[count];
        assert(column->size() == expected_ids.size());
        for (size_t i = 0; i < expected_ids.size(); i++) {
            const Entry* entry = column->at(i);
            assert(entry->get_read_id() == expected_ids[i]);
        }
        count++;
    }
    assert_msg(true, "ColumnIterator", "Column advancement test.");
    assert_msg(count == 2, "ColumnIterator", "2 columns iterated.");

    delete read_set;
    delete column_iterator;

}

#endif // TEST_COLUMNITERATOR_H