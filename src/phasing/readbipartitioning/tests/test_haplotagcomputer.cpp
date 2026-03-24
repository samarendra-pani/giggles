#include "test_haplotagcomputer.h"

void test_position_to_index() {

    ReadSet* superreads = mock_superreads();
    assert (superreads->size() == 2); // two superreads represeting the two haplotypes
    Read* superread0 = superreads->get(0);
    Read* superread1 = superreads->get(1);
    // Create a mapping from position to index for superreads
    std::unordered_map<uint32_t, uint32_t> position_to_index;
    for (uint32_t i = 0; i < superread0->getVariantCount(); ++i) {
        uint32_t pos = superread0->getPosition(i);
        assert(pos == superread1->getPosition(i));
        position_to_index[pos] = i;
    }
    for (uint32_t i = 0; i < superread0->getVariantCount(); ++i) {
        assert(position_to_index[(i+1)*100] == i);
    }
    assert_msg(true, "HaplotagComputer", "Position to index map.");
    delete superreads;
}

void test_haplotag_selected_reads() {
    std::vector<bool>* partitioning = new std::vector<bool>();
    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(true);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(true);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(false);
    
    partitioning->push_back(true);
    partitioning->push_back(false);
    partitioning->push_back(false);

    haplotag_selected_reads(read_set, partitioning);

    assert(read1->hasHaplotag());
    assert(read1->getHaplotag());
    assert(read2->hasHaplotag());
    assert(!read2->getHaplotag());
    assert(!read3->hasHaplotag());

    assert_msg(true, "HaplotagComputer", "Haplotags for selected reads.");

    delete read_set;
    delete partitioning;
}

void test_distance_calculation_from_superreads() {
    ReadSet* superreads = mock_superreads();
    assert (superreads->size() == 2); // two superreads represeting the two haplotypes
    Read* superread0 = superreads->get(0);
    Read* superread1 = superreads->get(1);
    // Create a mapping from position to index for superreads
    std::unordered_map<uint32_t, uint32_t> position_to_index;
    for (uint32_t i = 0; i < superread0->getVariantCount(); ++i) {
        uint32_t pos = superread0->getPosition(i);
        assert(pos == superread1->getPosition(i));
        position_to_index[pos] = i;
    }

    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(false);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(false);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(false);
    Read* read4 = new Read("read4", 60, 0); read_set->add(read4); read4->setSelected(false);
    Read* read5 = new Read("read5", 60, 0); read_set->add(read5); read5->setSelected(false);
    Read* read6 = new Read("read6", 60, 0); read_set->add(read6); read6->setSelected(false);
    Read* read7 = new Read("read7", 60, 0); read_set->add(read7); read7->setSelected(false);
    
    {
        read1->addVariant(100, std::vector<uint32_t>{}, Entry::ALLELE1);
        read1->addVariant(200, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(300, std::vector<uint32_t>{}, Entry::ALLELE2);
        
        read2->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read2->addVariant(400, std::vector<uint32_t>{}, Entry::ALLELE2);
        read2->addVariant(500, std::vector<uint32_t>{}, Entry::BLANK);
        read2->addVariant(600, std::vector<uint32_t>{}, Entry::ALLELE1);

        read3->addVariant(100, std::vector<uint32_t>{}, Entry::ALLELE2);
        read3->addVariant(200, std::vector<uint32_t>{}, Entry::ALLELE1);
        read3->addVariant(300, std::vector<uint32_t>{}, Entry::ALLELE2);
        read3->addVariant(400, std::vector<uint32_t>{}, Entry::ALLELE2);

        read4->addVariant(900, std::vector<uint32_t>{}, Entry::ALLELE2);
        read4->addVariant(1000, std::vector<uint32_t>{}, Entry::BLANK);
        read4->addVariant(1100, std::vector<uint32_t>{}, Entry::ALLELE1);
        read4->addVariant(1500, std::vector<uint32_t>{}, Entry::ALLELE1);

        read5->addVariant(1200, std::vector<uint32_t>{}, Entry::ALLELE1);
        read5->addVariant(1300, std::vector<uint32_t>{}, Entry::ALLELE2);
        read5->addVariant(1400, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read5->addVariant(1500, std::vector<uint32_t>{}, Entry::ALLELE1);
        
        read6->addVariant(1200, std::vector<uint32_t>{}, Entry::ALLELE1);
        read6->addVariant(1300, std::vector<uint32_t>{}, Entry::ALLELE1);
        read6->addVariant(1400, std::vector<uint32_t>{}, Entry::ALLELE2);
        read6->addVariant(1500, std::vector<uint32_t>{}, Entry::ALLELE2);

        read7->addVariant(1000, std::vector<uint32_t>{}, Entry::BLANK);
    }

    assert(calculate_distance_from_superread(read1, superread0, position_to_index) == 1);
    assert(calculate_distance_from_superread(read1, superread1, position_to_index) == 2);

    assert(calculate_distance_from_superread(read2, superread0, position_to_index) == 0);
    assert(calculate_distance_from_superread(read2, superread1, position_to_index) == 2);

    assert(calculate_distance_from_superread(read3, superread0, position_to_index) == 2);
    assert(calculate_distance_from_superread(read3, superread1, position_to_index) == 2);

    assert(calculate_distance_from_superread(read4, superread0, position_to_index) == 3);
    assert(calculate_distance_from_superread(read4, superread1, position_to_index) == 0);

    assert(calculate_distance_from_superread(read5, superread0, position_to_index) == 1);
    assert(calculate_distance_from_superread(read5, superread1, position_to_index) == 1);

    assert(calculate_distance_from_superread(read6, superread0, position_to_index) == 0);
    assert(calculate_distance_from_superread(read6, superread1, position_to_index) == 1);

    assert(calculate_distance_from_superread(read7, superread0, position_to_index) == 0);
    assert(calculate_distance_from_superread(read7, superread1, position_to_index) == 0);

    assert_msg(true, "HaplotagComputer", "Distance of read from superread.");
    
    delete superreads;
    delete read_set;
}

void test_haplotag_unselected_reads() {

    ReadSet* superreads = mock_superreads();
    assert (superreads->size() == 2); // two superreads represeting the two haplotypes
    Read* superread0 = superreads->get(0);
    Read* superread1 = superreads->get(1);
    std::vector<uint32_t> accessible_positions = {100, 200, 300, 400, 600, 700, 800, 900, 1100, 1200, 1300, 1400, 1500};
    // Create a mapping from position to index for superreads
    std::unordered_map<uint32_t, uint32_t> position_to_index;
    for (uint32_t i = 0; i < superread0->getVariantCount(); ++i) {
        uint32_t pos = superread0->getPosition(i);
        assert(pos == superread1->getPosition(i));
        position_to_index[pos] = i;
    }

    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(true);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(true);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(true);
    Read* read4 = new Read("read4", 60, 0); read_set->add(read4); read4->setSelected(true);
    Read* read5 = new Read("read5", 60, 0); read_set->add(read5); read5->setSelected(false);
    Read* read6 = new Read("read6", 60, 0); read_set->add(read6); read6->setSelected(false);
    Read* read7 = new Read("read7", 60, 0); read_set->add(read7); read7->setSelected(false);
    Read* read8 = new Read("read8", 60, 0); read_set->add(read8); read8->setSelected(false);
    Read* read9 = new Read("read9", 60, 0); read_set->add(read9); read9->setSelected(false);

    /**
     * First we create the phaseblocks using selected reads (reads 1, 2, 3, and 4)
     * 
     * phaseblocks are:
     * 100 -> {100, 400, 600}
     * 800 -> {800, 900, 1100}
     * 1500 -> {1500}
     */
    {
        read1->addVariant(100, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(200, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(400, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);

        read2->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read2->addVariant(400, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read2->addVariant(500, std::vector<uint32_t>{}, Entry::BLANK);
        read2->addVariant(600, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);

        read3->addVariant(800, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read3->addVariant(900, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read3->addVariant(1100, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);

        read4->addVariant(1100, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read4->addVariant(1200, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read4->addVariant(1300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);

        /** read5 should have no PS and no HP tag */
        read5->addVariant(600, std::vector<uint32_t>{}, Entry::ALLELE1);
        read5->addVariant(700, std::vector<uint32_t>{}, Entry::ALLELE2);
        read5->addVariant(800, std::vector<uint32_t>{}, Entry::ALLELE2);

        /** read6 should have PS=800 and HP=0 */
        read6->addVariant(900, std::vector<uint32_t>{}, Entry::ALLELE1);
        read6->addVariant(1000, std::vector<uint32_t>{}, Entry::BLANK);
        read6->addVariant(1100, std::vector<uint32_t>{}, Entry::ALLELE2);
        read6->addVariant(1200, std::vector<uint32_t>{}, Entry::ALLELE1);
        read6->addVariant(1300, std::vector<uint32_t>{}, Entry::ALLELE1);

        /** read7 should have no PS and no HP tag */
        read7->addVariant(1000, std::vector<uint32_t>{}, Entry::BLANK);

        /** read8 is equidistant from superread 0 and 1. PS=1500 and no HP. */
        read8->addVariant(1300, std::vector<uint32_t>{}, Entry::ALLELE2);
        read8->addVariant(1400, std::vector<uint32_t>{}, Entry::ALLELE2);
        read8->addVariant(1500, std::vector<uint32_t>{}, Entry::ALLELE1);

        /** read9 should have PS=1500 and HP=0 */
        read9->addVariant(1300, std::vector<uint32_t>{}, Entry::ALLELE2);
        read9->addVariant(1400, std::vector<uint32_t>{}, Entry::ALLELE2);
        read9->addVariant(1500, std::vector<uint32_t>{}, Entry::ALLELE2);
    }

    compute_phasesets(&accessible_positions, read_set, superreads);
    /** checking phasesets. */
    {
        assert(read1->hasPhaseSet());
        assert(read1->getPhaseSet() == 100);
        assert(read2->hasPhaseSet());
        assert(read2->getPhaseSet() == 100);
        assert(read3->hasPhaseSet());
        assert(read3->getPhaseSet() == 800);
        assert(read4->hasPhaseSet());
        assert(read4->getPhaseSet() == 800);
        assert(!read5->hasPhaseSet());
        assert(read6->hasPhaseSet());
        assert(read6->getPhaseSet() == 800);
        assert(!read7->hasPhaseSet());
        assert(read8->hasPhaseSet());
        assert(read8->getPhaseSet() == 1500);
        assert(read9->hasPhaseSet());
        assert(read9->getPhaseSet() == 1500);
    }
    
    std::vector<bool>* partitioning = new std::vector<bool>();
    /** adding data for partitioning. */
    {
        partitioning->push_back(true); // read1 has HP=1
        partitioning->push_back(false); // read2 has HP=0
        partitioning->push_back(false); // read3 has HP=0
        partitioning->push_back(true); // read4 has HP=1
        partitioning->push_back(false); //read5 is not selected. So false by default
        partitioning->push_back(false); //read6 is not selected. So false by default
        partitioning->push_back(false); //read7 is not selected. So false by default
        partitioning->push_back(false); //read8 is not selected. So false by default
        partitioning->push_back(false); //read9 is not selected. So false by default
    }
    haplotag_selected_reads(read_set, partitioning);
    /** checking haplotags before tagging unselected reads. */
    {
        assert(read1->hasHaplotag());
        assert(read1->getHaplotag());
        assert(read2->hasHaplotag());
        assert(!read2->getHaplotag());
        assert(read3->hasHaplotag());
        assert(!read3->getHaplotag());
        assert(read4->hasHaplotag());
        assert(read4->getHaplotag());
        assert(!read5->hasHaplotag());
        assert(!read6->hasHaplotag()); 
        assert(!read7->hasHaplotag()); 
        assert(!read8->hasHaplotag());
        assert(!read9->hasHaplotag());
    }

    haplotag_unselected_reads(read_set, superreads);
    {
        assert(read1->hasHaplotag());
        assert(read1->getHaplotag());
        assert(read2->hasHaplotag());
        assert(!read2->getHaplotag());
        assert(read3->hasHaplotag());
        assert(!read3->getHaplotag());
        assert(read4->hasHaplotag());
        assert(read4->getHaplotag());
        assert(!read5->hasHaplotag());
        assert(read6->hasHaplotag()); 
        assert(!read6->getHaplotag());
        assert(!read7->hasHaplotag()); 
        assert(!read8->hasHaplotag());
        assert(read9->hasHaplotag());
        assert(!read9->getHaplotag());  
    }
    assert_msg(true, "HaplotagComputer", "Haplotags for unselected reads.");
    
    delete superreads;
    delete read_set;
    delete partitioning;
}

void test_haplotagcomputer() {

    test_position_to_index();
    test_haplotag_selected_reads();
    test_distance_calculation_from_superreads();
    test_haplotag_unselected_reads();   
}