#include "test_readset.h"

void test_pos_to_entry_map() {
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_2();

    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(false);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(true);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(true);
    Read* read4 = new Read("read4", 60, 0); read_set->add(read4); read4->setSelected(true);
    Read* read5 = new Read("read5", 60, 0); read_set->add(read5); read5->setSelected(true);
    Read* read6 = new Read("read6", 60, 0); read_set->add(read6); read6->setSelected(false);

    {   
        read1->addVariant(100, std::vector<uint32_t>{10, 20});

        read2->addVariant(200, std::vector<uint32_t>{20, 5, 10});
        
        read3->addVariant(200, std::vector<uint32_t>{15, 8, 10});
        
        read4->addVariant(400, std::vector<uint32_t>{20, 5, 2, 1});

        read5->addVariant(400, std::vector<uint32_t>{10, 15, 1, 2});
        read5->addVariant(500, std::vector<uint32_t>{20, 5});

        read6->addVariant(500, std::vector<uint32_t>{15, 8});
    }
    read_set->initialize();

    assert(read_set->TEST_get_pos_to_entry_map(100).size() == 1);
    assert(read_set->TEST_get_pos_to_entry_map(100)[0]->get_read_id() == 0);
    assert(read_set->TEST_get_pos_to_entry_map(200).size() == 2);
    assert(read_set->TEST_get_pos_to_entry_map(200)[0]->get_read_id() == 1);
    assert(read_set->TEST_get_pos_to_entry_map(200)[1]->get_read_id() == 2);
    assert(read_set->TEST_get_pos_to_entry_map(300).size() == 0);
    assert(read_set->TEST_get_pos_to_entry_map(400).size() == 2);
    assert(read_set->TEST_get_pos_to_entry_map(400)[0]->get_read_id() == 3);
    assert(read_set->TEST_get_pos_to_entry_map(400)[1]->get_read_id() == 4);
    assert(read_set->TEST_get_pos_to_entry_map(500).size() == 2);
    assert(read_set->TEST_get_pos_to_entry_map(500)[0]->get_read_id() == 4);
    assert(read_set->TEST_get_pos_to_entry_map(500)[1]->get_read_id() == 5);
    
    assert_msg(true, "ReadSet", "Position to Entry Map.");
}

void test_readset_sorting() {
    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(false);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(true);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(true);
    Read* read4 = new Read("read4", 60, 0); read_set->add(read4); read4->setSelected(true);
    Read* read5 = new Read("read5", 60, 0); read_set->add(read5); read5->setSelected(true);
    Read* read6 = new Read("read6", 60, 0); read_set->add(read6); read6->setSelected(false);

    {   
        read1->addVariant(100, std::vector<uint32_t>{10, 20});

        read3->addVariant(200, std::vector<uint32_t>{15, 8, 10});
        
        read4->addVariant(400, std::vector<uint32_t>{20, 5, 2, 1});

        read5->addVariant(300, std::vector<uint32_t>{10, 15, 1, 2});
        read5->addVariant(500, std::vector<uint32_t>{20, 5});

        read6->addVariant(200, std::vector<uint32_t>{15, 8});
    }
    read_set->initialize();
    assert(read2->getID() == 0);
    assert(read1->getID() == 1);
    assert(read3->getID() == 2);
    assert(read6->getID() == 3);
    assert(read5->getID() == 4);
    assert(read4->getID() == 5);

    assert_msg(true, "ReadSet", "Sorting of Reads.");
}

void test_readset() {
    test_pos_to_entry_map();
    test_readset_sorting();
}