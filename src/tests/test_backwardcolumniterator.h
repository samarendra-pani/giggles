#ifndef TEST_BACKWARDCOLUMNITERATOR_H
#define TEST_BACKWARDCOLUMNITERATOR_H

#include "../backwardcolumniterator.h"
#include "../tests_data.h"
#include <cassert>


void test_backwardcolumniterator() {
    
    std::string prefix = "BackwardColumnIterator";
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_1();
    ReadSet* read_set = mock_readset_1();

    BackwardColumnIterator* column_iterator = new BackwardColumnIterator(*read_set, &variant_info_table);
    
    assert_msg(column_iterator->get_column_count() == 2, prefix, "Column count should be 2.");
    assert_msg(column_iterator->get_read_count() == read_set->size(), prefix, "Read count should match read set size.");
    const std::vector<uint32_t>* positions = column_iterator->get_positions();
    assert_msg(positions->size() == 2, prefix, "Positions size should be 2.");
    assert(positions->at(0) == 100);
    assert(positions->at(1) == 200);

    std::cout << "[BackwardColumnIterator] Testing advancement..." << std::endl;
    
     std::vector<std::vector<uint32_t>> expected_entry_id = {
        {12, 14, 16, 17, 19, 20, 21, 22, 24, 25, 26, 27, 28, 29},
        {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27}
    };

    uint32_t count = 0;
    while (column_iterator->has_next()) {
        std::unique_ptr<std::vector<const Entry*> > column = column_iterator->get_next();
        std::cout << "[BackwardColumnIterator] Testing column " << count << "..." << std::endl;
        // check expected entries
        const std::vector<uint32_t>& expected_ids = expected_entry_id[count];
        assert(column->size() == expected_ids.size());
        for (size_t i = 0; i < expected_ids.size(); i++) {
            const Entry* entry = column->at(i);
            assert(entry->get_read_id() == expected_ids[i]);
        }
        count++;
        assert_msg(true, prefix, "Advanced to next column successfully.");
    }
    assert_msg(count == 2, prefix, "Total columns iterated should be 2.");

    delete read_set;
    delete column_iterator;

}

#endif // TEST_BACKWARDCOLUMNITERATOR_H