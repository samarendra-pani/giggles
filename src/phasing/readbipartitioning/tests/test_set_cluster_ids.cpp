#include "test_set_cluster_ids.h"

void test_singleposition_blankentries() {
    ReadSet* superreads = mock_superreads();
    Read* superread0 = superreads->get(0);
    Read* superread1 = superreads->get(1);
    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(false);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(false);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(false);
    Read* read4 = new Read("read4", 60, 0); read_set->add(read4); read4->setSelected(false);
    Read* read5 = new Read("read5", 60, 0); read_set->add(read5); read5->setSelected(false);
    Read* read6 = new Read("read6", 60, 0); read_set->add(read6); read6->setSelected(false);
    Read* read7 = new Read("read7", 60, 0); read_set->add(read7); read7->setSelected(false);
    Read* read8 = new Read("read8", 60, 0); read_set->add(read8); read8->setSelected(false);
    Read* read9 = new Read("read9", 60, 0); read_set->add(read9); read9->setSelected(false);
    Read* read10 = new Read("read10", 60, 0); read_set->add(read10); read10->setSelected(false);
    std::vector<uint32_t> accessible_positions = {100, 200, 300, 400, 600, 700, 800, 900, 1100, 1200, 1300, 1400, 1500};
    // Create a mapping from position to index for superreads
    std::unordered_map<uint32_t, uint32_t> position_to_index;
    for (uint32_t i = 0; i < superread0->getVariantCount(); ++i) {
        uint32_t pos = superread0->getPosition(i);
        assert(pos == superread1->getPosition(i));
        position_to_index[pos] = i;
    }
    

    read1->addVariant(500, std::vector<uint32_t>{});
    read2->addVariant(500, std::vector<uint32_t>{});
    read3->addVariant(500, std::vector<uint32_t>{});
    read4->addVariant(500, std::vector<uint32_t>{});
    read5->addVariant(500, std::vector<uint32_t>{});
    read6->addVariant(500, std::vector<uint32_t>{});
    read7->addVariant(500, std::vector<uint32_t>{});
    read8->addVariant(500, std::vector<uint32_t>{});
    read9->addVariant(500, std::vector<uint32_t>{});
    read10->addVariant(500, std::vector<uint32_t>{});

    read1->setID(0);
    read2->setID(1);
    read3->setID(2);
    read4->setID(3);
    read5->setID(4);
    read6->setID(5);
    read7->setID(6);
    read8->setID(7);
    read9->setID(8);
    read10->setID(9);

    compute_phasesets(&accessible_positions, read_set, superreads);
    std::vector<bool>* partitioning = new std::vector<bool>();
    for (uint32_t i = 0; i < read_set->size(); i++) {
        partitioning->push_back(false);
    }
    haplotag_selected_reads(read_set, partitioning);
    haplotag_unselected_reads(read_set, superreads);
    set_read_cluster_ids(read_set);
    
    for (uint32_t i = 0; i < read_set->size(); i++) {
        assert(!read_set->get(i)->hasHaplotag());
        assert(!read_set->get(i)->hasPhaseSet());
        assert(!read_set->get(i)->isClustered());
        assert(read_set->get(i)->getClusterID() == i);
        assert(read_set->get(i)->getConstrainedClusterID() == (uint32_t)-1);
    }

    assert_msg(true, "SetClusterIDs", "Unselected single-position BLANK reads mapped to BLANK position.");
}

void test_multiposition_entries() {
    ReadSet* superreads = mock_superreads();
    Read* superread0 = superreads->get(0);
    Read* superread1 = superreads->get(1);
    ReadSet* read_set = new ReadSet();

    /**
     * Set of selected reads determine the phase blocks.
     * Haplotag determined by optimal_partitioning
     * 
     * Set of unselected reads are fit into phase blocks.
     * Haplotags determined by distance to superreads.
     */
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(true);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(false);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(false);
    Read* read4 = new Read("read4", 60, 0); read_set->add(read4); read4->setSelected(false);
    Read* read5 = new Read("read5", 60, 0); read_set->add(read5); read5->setSelected(false);
    Read* read6 = new Read("read6", 60, 0); read_set->add(read6); read6->setSelected(false);
    Read* read7 = new Read("read7", 60, 0); read_set->add(read7); read7->setSelected(true);
    Read* read8 = new Read("read8", 60, 0); read_set->add(read8); read8->setSelected(false);
    Read* read9 = new Read("read9", 60, 0); read_set->add(read9); read9->setSelected(false);
    Read* read10 = new Read("read10", 60, 0); read_set->add(read10); read10->setSelected(true);
    Read* read11 = new Read("read11", 60, 0); read_set->add(read11); read11->setSelected(false);
    Read* read12 = new Read("read12", 60, 0); read_set->add(read12); read12->setSelected(false);
    Read* read13 = new Read("read13", 60, 0); read_set->add(read13); read13->setSelected(true);
    

    {
        /* Selected Reads */
        /* PS: 100, HP: 0 */
        read1->addVariant(100, std::vector<uint32_t>{}, Entry::ALLELE1);
        read1->addVariant(200, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(400, std::vector<uint32_t>{}, Entry::ALLELE2);
        /* PS: 600, HP: 1 */
        read7->addVariant(600, std::vector<uint32_t>{}, Entry::ALLELE2); 
        read7->addVariant(700, std::vector<uint32_t>{}, Entry::ALLELE2);
        read7->addVariant(800, std::vector<uint32_t>{}, Entry::ALLELE1);
        read7->addVariant(900, std::vector<uint32_t>{}, Entry::ALLELE1);
        /* PS: 1100, HP: 0 */
        read10->addVariant(1000, std::vector<uint32_t>{}, Entry::BLANK);
        read10->addVariant(1100, std::vector<uint32_t>{}, Entry::ALLELE2);
        read10->addVariant(1200, std::vector<uint32_t>{}, Entry::ALLELE1);
        read10->addVariant(1300, std::vector<uint32_t>{}, Entry::ALLELE1);
        /* PS: 1500, HP: 0 */
        read13->addVariant(1300, std::vector<uint32_t>{}, Entry::ALLELE1);
        read13->addVariant(1400, std::vector<uint32_t>{}, Entry::ALLELE2);
        read13->addVariant(1500, std::vector<uint32_t>{}, Entry::ALLELE2);

        /**
         * Resulting phasesets are
         * 100 -> [100, 400]
         * 600 -> [600, 800, 900]
         * 1100 -> [1100]
         * 1500 -> [1500]
         */
    }

    {
        /* Unselected Reads */
        /* PS: 100, HP: 0 */
        read2->addVariant(100, std::vector<uint32_t>{}, Entry::ALLELE1);
        read2->addVariant(200, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read2->addVariant(300, std::vector<uint32_t>{}, Entry::ALLELE1);
        /* PS: 100, HP: 1 */
        read3->addVariant(200, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read3->addVariant(300, std::vector<uint32_t>{}, Entry::ALLELE1);
        read3->addVariant(400, std::vector<uint32_t>{}, Entry::ALLELE1);
        /* PS: ??, HP: ?? */
        read4->addVariant(300, std::vector<uint32_t>{}, Entry::ALLELE1);
        /* PS: ??, HP: ?? */
        read5->addVariant(400, std::vector<uint32_t>{}, Entry::ALLELE1);
        read5->addVariant(500, std::vector<uint32_t>{}, Entry::BLANK);
        read5->addVariant(600, std::vector<uint32_t>{}, Entry::ALLELE2);
        /* PS: 600, HP: 0 */
        read6->addVariant(500, std::vector<uint32_t>{}, Entry::BLANK);
        read6->addVariant(600, std::vector<uint32_t>{}, Entry::ALLELE1);
        read6->addVariant(700, std::vector<uint32_t>{}, Entry::ALLELE1);
        read6->addVariant(800, std::vector<uint32_t>{}, Entry::ALLELE2);
        /* PS: 600, HP: 1 */
        read8->addVariant(600, std::vector<uint32_t>{}, Entry::ALLELE2);
        read8->addVariant(700, std::vector<uint32_t>{}, Entry::ALLELE2);
        read8->addVariant(800, std::vector<uint32_t>{}, Entry::ALLELE1);
        /* PS: 600, HP: 1 */
        read9->addVariant(800, std::vector<uint32_t>{}, Entry::ALLELE1);
        read9->addVariant(900, std::vector<uint32_t>{}, Entry::ALLELE2);
        read9->addVariant(1000, std::vector<uint32_t>{}, Entry::BLANK);
        /* PS: 1100. HP: ?? */
        read11->addVariant(1000, std::vector<uint32_t>{}, Entry::BLANK);
        read11->addVariant(1100, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        /* PS: 1100. HP: ?? */
        read12->addVariant(1100, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
    }
    
    std::vector<bool>* partitioning = new std::vector<bool>(read_set->size(), false);
    {
        partitioning->at(0) = false;   // for read1
        partitioning->at(6) = true;    // for read7
        partitioning->at(9) = false;   // for read10
        partitioning->at(12) = false;   // for read13
    }
    
    {
        /* Manually setting IDs to avoid hash function shennanigans. */
        read1->setID(0);
        read2->setID(1);
        read3->setID(2);
        read4->setID(3);
        read5->setID(4);
        read6->setID(5);
        read7->setID(6);
        read8->setID(7);
        read9->setID(8);
        read10->setID(9);
        read11->setID(10);
        read12->setID(11);
        read13->setID(12);
    }

    std::vector<uint32_t> accessible_positions = {100, 200, 300, 400, 600, 700, 800, 900, 1100, 1200, 1300, 1400, 1500};
    // Create a mapping from position to index for superreads
    std::unordered_map<uint32_t, uint32_t> position_to_index;
    for (uint32_t i = 0; i < superread0->getVariantCount(); ++i) {
        uint32_t pos = superread0->getPosition(i);
        assert(pos == superread1->getPosition(i));
        position_to_index[pos] = i;
    }
    
    compute_phasesets(&accessible_positions, read_set, superreads);
    haplotag_selected_reads(read_set, partitioning);
    haplotag_unselected_reads(read_set, superreads);
    set_read_cluster_ids(read_set);
    
    std::vector<bool> expected_has_hp = {true, true, true, false, false, true, true, true, true, true, false, false, true};
    std::vector<bool> expected_has_ps = {true, true, true, false, false, true, true, true, true, true, true, true, true};
    std::vector<bool> expected_is_clustered = {true, true, true, false, false, true, true, true, true, true, false, false, true};
    std::vector<bool> expected_is_constrained = {true, true, true, false, false, true, true, true, true, false, false, false, false};
    std::vector<uint32_t> expected_ps = {100, 100, 100, 0, 0, 600, 600, 600, 600, 1100, 1100, 1100, 1500};
    std::vector<bool> expected_hp = {false, false, true, false, false, false, true, true, true, false, false, false, false};
    std::vector<uint32_t> expected_cluster_id = {0, 0, 2, 3, 4, 5, 6, 6, 6, 9, 10, 11 , 12};
    std::vector<uint32_t> expected_constrained_id = {2, 2, 0, (uint32_t)-1, (uint32_t)-1, 6, 5, 5, 5, (uint32_t)-1, (uint32_t)-1, (uint32_t)-1, (uint32_t)-1};
    
    

    for (uint32_t i = 0; i < read_set->size(); i++) {
        /*std::cout << "ID: " << read_set->get(i)->getID()  << std::endl;
        std::cout << "Has HP: " << read_set->get(i)->hasHaplotag() << std::endl;
        if (read_set->get(i)->hasHaplotag()) {
            std::cout << "HP: " << read_set->get(i)->getHaplotag() << std::endl;
        }
        std::cout << "Has PS: " << read_set->get(i)->hasPhaseSet() << std::endl;
        if (read_set->get(i)->hasPhaseSet()) {
            std::cout << "PS: " << read_set->get(i)->getPhaseSet() << std::endl;
        }
        std::cout << "Is CS: " << read_set->get(i)->isClustered() << std::endl;
        std::cout << "Cluster ID: " << read_set->get(i)->getClusterID() << std::endl;
        std::cout << "Has Constrained ID: " << read_set->get(i)->hasConstrainedCluster() << std::endl;
        std::cout << "Constrained ID: " << read_set->get(i)->getConstrainedClusterID() << std::endl;
        std::cout << std::endl;*/
        
        assert(read_set->get(i)->hasHaplotag() == expected_has_hp[i]);
        assert(read_set->get(i)->hasPhaseSet() == expected_has_ps[i]);
        if (read_set->get(i)->hasHaplotag()) { assert(read_set->get(i)->getHaplotag() == expected_hp[i]); }
        if (read_set->get(i)->hasPhaseSet()) { assert(read_set->get(i)->getPhaseSet() == expected_ps[i]); }
        assert(read_set->get(i)->isClustered() == expected_is_clustered[i]);
        assert(read_set->get(i)->getClusterID() == expected_cluster_id[i]);
        assert(read_set->get(i)->hasConstrainedCluster() == expected_is_constrained[i]);
        assert(read_set->get(i)->getConstrainedClusterID() == expected_constrained_id[i]);
    }

    assert_msg(true, "SetClusterIDs", "Setting cluster IDs with mix of selected and unselected reads.");
}

void test_set_cluster_ids() {

    test_singleposition_blankentries();
    test_multiposition_entries();
}