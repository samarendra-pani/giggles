#include "test_phasingdptable.h"

void test_singleposition_allreadsselected() {
    std::vector<variant_information_t> variant_info_table;
    std::vector<uint32_t> position = {100};
    uint32_t ploidy = 2;
    std::vector<uint32_t> n_alleles = {2};
    std::vector<std::vector<int>> allele_references = {{1, 0, 1, 0}};
    std::vector<bool> is_sv_position = {false};
    variant_info_table.push_back(variant_information_t(position[0], ploidy, n_alleles[0], allele_references[0], is_sv_position[0]));

    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(true);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(true);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(true);
    Read* read4 = new Read("read4", 60, 0); read_set->add(read4); read4->setSelected(true);
    Read* read5 = new Read("read5", 60, 0); read_set->add(read5); read5->setSelected(true);
    Read* read6 = new Read("read6", 60, 0); read_set->add(read6); read6->setSelected(true);
    Read* read7 = new Read("read7", 60, 0); read_set->add(read7); read7->setSelected(true);
    Read* read8 = new Read("read8", 60, 0); read_set->add(read8); read8->setSelected(true);
    Read* read9 = new Read("read9", 60, 0); read_set->add(read9); read9->setSelected(true);
    Read* read10 = new Read("read10", 60, 0); read_set->add(read10); read10->setSelected(true);

    read1->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1
    read2->addVariant(100, std::vector<uint32_t>{10, 5}); // ALLELE2
    read3->addVariant(100, std::vector<uint32_t>{10, 10}); // EQUAL_SCORES
    read4->addVariant(100, std::vector<uint32_t>{10, 10}); // EQUAL_SCORES
    read5->addVariant(100, std::vector<uint32_t>{10, 5}); // ALLELE2
    read6->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1
    read7->addVariant(100, std::vector<uint32_t>{20, 10}); // ALLELE2
    read8->addVariant(100, std::vector<uint32_t>{10, 20}); // ALLELE1
    read9->addVariant(100, std::vector<uint32_t>{5, 10}); // ALLELE1
    read10->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1

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

    PhasingDPTable table = PhasingDPTable(read_set, &variant_info_table, true);
    
    std::vector<uint32_t> expected_cluster_ids = {0, 1, 0, 0, 1, 0, 1, 0, 0, 0};
    std::vector<bool> expected_cluster_status = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    std::vector<uint32_t> expected_constrained_cluster_ids = {1, 0, 1, 1, 0, 1, 0, 1, 1, 1};
    std::vector<bool> expected_constrained_cluster_status = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    assert(expected_cluster_ids.size() == read_set->size());
    for (uint32_t i = 0; i < read_set->size(); i++) {
        //std::cout << read_set->get(i)->getName() << "\t\t" << read_set->get(i)->isClustered() << "\t\t" << read_set->get(i)->getClusterID() <<"\t\t" << read_set->get(i)->getConstrainedClusterID() << std::endl;
        assert(expected_cluster_status[i] == read_set->get(i)->isClustered());
        assert(expected_cluster_ids[i] == read_set->get(i)->getClusterID());
        assert(expected_constrained_cluster_status[i] == read_set->get(i)->hasConstrainedCluster());
        if (expected_constrained_cluster_status[i]) {
            assert(expected_constrained_cluster_ids[i] == read_set->get(i)->getConstrainedClusterID());
        }
    }
    assert_msg(true, "PhasingDPTable", "Cluster IDs for single position variant table with all reads selected.");
    delete read_set;
}

void test_singleposition_somereadsselected() {
    std::vector<variant_information_t> variant_info_table;
    std::vector<uint32_t> position = {100};
    uint32_t ploidy = 2;
    std::vector<uint32_t> n_alleles = {2};
    std::vector<std::vector<int>> allele_references = {{1, 0, 1, 0}};
    std::vector<bool> is_sv_position = {false};
    variant_info_table.push_back(variant_information_t(position[0], ploidy, n_alleles[0], allele_references[0], is_sv_position[0]));

    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(true);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(true);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(false);
    Read* read4 = new Read("read4", 60, 0); read_set->add(read4); read4->setSelected(false);
    Read* read5 = new Read("read5", 60, 0); read_set->add(read5); read5->setSelected(false);
    Read* read6 = new Read("read6", 60, 0); read_set->add(read6); read6->setSelected(true);
    Read* read7 = new Read("read7", 60, 0); read_set->add(read7); read7->setSelected(true);
    Read* read8 = new Read("read8", 60, 0); read_set->add(read8); read8->setSelected(true);
    Read* read9 = new Read("read9", 60, 0); read_set->add(read9); read9->setSelected(false);
    Read* read10 = new Read("read10", 60, 0); read_set->add(read10); read10->setSelected(true);

    read1->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1
    read2->addVariant(100, std::vector<uint32_t>{10, 5}); // ALLELE2
    read3->addVariant(100, std::vector<uint32_t>{10, 10}); // EQUAL_SCORES
    read4->addVariant(100, std::vector<uint32_t>{10, 10}); // EQUAL_SCORES
    read5->addVariant(100, std::vector<uint32_t>{10, 5}); // ALLELE2
    read6->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1
    read7->addVariant(100, std::vector<uint32_t>{20, 10}); // ALLELE2
    read8->addVariant(100, std::vector<uint32_t>{10, 20}); // ALLELE1
    read9->addVariant(100, std::vector<uint32_t>{5, 10}); // ALLELE1
    read10->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1

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

    PhasingDPTable table = PhasingDPTable(read_set, &variant_info_table, true);
    
    std::vector<uint32_t> expected_cluster_ids = {0, 1, 2, 3, 1, 0, 1, 0, 0, 0};
    std::vector<bool> expected_cluster_status = {1, 1, 0, 0, 1, 1, 1, 1, 1, 1};
    std::vector<uint32_t> expected_constrained_cluster_ids = {1, 0, -1u, -1u, 0, 1, 0, 1, 1, 1};
    std::vector<bool> expected_constrained_cluster_status = {1, 1, 0, 0, 1, 1, 1, 1, 1, 1};
    assert(expected_cluster_ids.size() == read_set->size());
    for (uint32_t i = 0; i < read_set->size(); i++) {
        //std::cout << read_set->get(i)->getName() << "\t\t" << read_set->get(i)->isClustered() << "\t\t" << read_set->get(i)->getClusterID() <<"\t\t" << read_set->get(i)->getConstrainedClusterID() << std::endl;
        assert(expected_cluster_status[i] == read_set->get(i)->isClustered());
        assert(expected_cluster_ids[i] == read_set->get(i)->getClusterID());
        assert(expected_constrained_cluster_status[i] == read_set->get(i)->hasConstrainedCluster());
        if (expected_constrained_cluster_status[i]) {
            assert(expected_constrained_cluster_ids[i] == read_set->get(i)->getConstrainedClusterID());
        }
    }
    assert_msg(true, "PhasingDPTable", "Cluster IDs for single position variant table with some reads selected.");
    delete read_set;
}

void test_singleposition_sv() {
    std::vector<variant_information_t> variant_info_table;
    std::vector<uint32_t> position = {100};
    uint32_t ploidy = 2;
    std::vector<uint32_t> n_alleles = {3};
    std::vector<std::vector<int>> allele_references = {{1, 0, 1, 2}};
    std::vector<bool> is_sv_position = {true};
    variant_info_table.push_back(variant_information_t(position[0], ploidy, n_alleles[0], allele_references[0], is_sv_position[0]));

    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(true);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(true);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(true);
    Read* read4 = new Read("read4", 60, 0); read_set->add(read4); read4->setSelected(true);
    Read* read5 = new Read("read5", 60, 0); read_set->add(read5); read5->setSelected(true);
    Read* read6 = new Read("read6", 60, 0); read_set->add(read6); read6->setSelected(true);
    Read* read7 = new Read("read7", 60, 0); read_set->add(read7); read7->setSelected(true);
    Read* read8 = new Read("read8", 60, 0); read_set->add(read8); read8->setSelected(true);
    Read* read9 = new Read("read9", 60, 0); read_set->add(read9); read9->setSelected(true);
    Read* read10 = new Read("read10", 60, 0); read_set->add(read10); read10->setSelected(true);

    /** In the end, 0 and 1 will be the active position. */
    read1->addVariant(100, std::vector<uint32_t>{10, 90, 50}); // ALLELE1
    read2->addVariant(100, std::vector<uint32_t>{10, 5, 50}); // ALLELE2
    read3->addVariant(100, std::vector<uint32_t>{10, 10, 50}); // EQUAL_SCORES
    read4->addVariant(100, std::vector<uint32_t>{10, 10, 50}); // EQUAL_SCORES
    read5->addVariant(100, std::vector<uint32_t>{10, 5, 50}); // ALLELE2
    read6->addVariant(100, std::vector<uint32_t>{10, 90, 50}); // ALLELE1
    read7->addVariant(100, std::vector<uint32_t>{20, 10, 50}); // ALLELE2
    read8->addVariant(100, std::vector<uint32_t>{10, 20, 50}); // ALLELE1
    read9->addVariant(100, std::vector<uint32_t>{5, 10, 50}); // ALLELE1
    read10->addVariant(100, std::vector<uint32_t>{10, 90, 50}); // ALLELE1

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

    /** First round of phasing should not consider SVs */
    PhasingDPTable* table = new PhasingDPTable(read_set, &variant_info_table, true);
    
    std::vector<uint32_t> expected_cluster_ids = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<bool> expected_cluster_status = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<uint32_t> expected_constrained_cluster_ids = {-1u, -1u, -1u, -1u, -1u, -1u, -1u, -1u, -1u, -1u};
    std::vector<bool> expected_constrained_cluster_status = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    assert(expected_cluster_ids.size() == read_set->size());
    for (uint32_t i = 0; i < read_set->size(); i++) {
        assert(expected_cluster_status[i] == read_set->get(i)->isClustered());
        assert(expected_cluster_ids[i] == read_set->get(i)->getClusterID());
        assert(expected_constrained_cluster_status[i] == read_set->get(i)->hasConstrainedCluster());
        if (expected_constrained_cluster_status[i]) {
            assert(expected_constrained_cluster_ids[i] == read_set->get(i)->getConstrainedClusterID());
        }
    }

    delete table;
    read_set->resetTags();
    variant_info_table[0].genotype_likelihoods.increment_by_index(0, 0.5L);  // selected
    variant_info_table[0].genotype_likelihoods.increment_by_index(1, 0.3L);  // selected
    variant_info_table[0].genotype_likelihoods.increment_by_index(2, 0.15L); // selected
    variant_info_table[0].genotype_likelihoods.increment_by_index(3, 0.05L);
    std::vector<uint32_t> selected_genotype_indices = variant_info_table[0].genotype_likelihoods.select_genotypes();
	variant_info_table[0].update_active_alleles(2, selected_genotype_indices);
    
    /** Subsequent rounds of phasing should consider SVs if they are bi-allelic */
    table = new PhasingDPTable(read_set, &variant_info_table, false);
    
    expected_cluster_ids = {0, 1, 0, 0, 1, 0, 1, 0, 0, 0};
    expected_cluster_status = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    expected_constrained_cluster_ids = {1, 0, 1, 1, 0, 1, 0, 1, 1, 1};
    expected_constrained_cluster_status = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};assert(expected_cluster_ids.size() == read_set->size());
    for (uint32_t i = 0; i < read_set->size(); i++) {
        //std::cout << read_set->get(i)->getName() << "\t\t" << read_set->get(i)->isClustered() << "\t\t" << read_set->get(i)->getClusterID() <<"\t\t" << read_set->get(i)->getConstrainedClusterID() << std::endl;
        assert(expected_cluster_status[i] == read_set->get(i)->isClustered());
        assert(expected_cluster_ids[i] == read_set->get(i)->getClusterID());
        assert(expected_constrained_cluster_status[i] == read_set->get(i)->hasConstrainedCluster());
        if (expected_constrained_cluster_status[i]) {
            assert(expected_constrained_cluster_ids[i] == read_set->get(i)->getConstrainedClusterID());
        }
    }
    delete table;
    delete read_set;
    assert_msg(true, "PhasingDPTable", "Cluster IDs for single SV position variant table.");
}

void test_multiposition() {
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_2();
    /**
     * Variant Information Table:
     * ........................| variant 1 | variant 2 | variant 3 | variant 4 | variant 5
     * ------------------------|-----------|-----------|-----------|-----------|-----------
     * Position                | 100       | 200       | 300       | 400       | 500
     * Number of Alleles       | 2         | 3         | 2         | 4         | 2
     * Is a SV?                | No        | Yes       | No        | Yes       | No
     * -------------------------------------------------------------------------------------
     * Haplotype 1 Allele      | 1         | 0         | 1         | 0         | 0
     * Haplotype 2 Allele      | 0         | 1         | 1         | 1         | 0
     * Haplotype 3 Allele      | 1         | 0         | 0         | 3         | 1
     * Haplotype 4 Allele      | 0         | 2         | 0         | 2         | 0
     */

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
    
    /* First round of phasing */
    {
        PhasingDPTable table = PhasingDPTable(read_set, &variant_info_table, true);

        std::vector<uint32_t> expected_cluster_ids = {0, 1, 2, 3, 4, 5, 6, 7, 6, 7, 10, 11, 11};
        std::vector<bool> expected_cluster_status = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1};
        std::vector<uint32_t> expected_constrained_cluster_ids = {-1u, -1u, -1u, -1u, -1u, -1u, 7, 6, 7, 6, 11, 10, 10};
        std::vector<bool> expected_constrained_cluster_status = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1};
        assert(expected_cluster_ids.size() == read_set->size());
        for (uint32_t i = 0; i < read_set->size(); i++) {
            //std::cout << read_set->get(i)->getName() << "\t\t" << read_set->get(i)->isClustered() << "\t\t" << read_set->get(i)->getClusterID() <<"\t\t" << read_set->get(i)->getConstrainedClusterID() << std::endl;
            assert(expected_cluster_status[i] == read_set->get(i)->isClustered());
            assert(expected_cluster_ids[i] == read_set->get(i)->getClusterID());
            assert(expected_constrained_cluster_status[i] == read_set->get(i)->hasConstrainedCluster());
            if (expected_constrained_cluster_status[i]) {
                assert(expected_constrained_cluster_ids[i] == read_set->get(i)->getConstrainedClusterID());
            }
        }
    }

    /* Updating GenotypeLikelihoods for variant at 200 so that only allele 1 and 2 are chosen*/
    variant_info_table[1].genotype_likelihoods.increment_by_index(2, 0.5L); // 2 -> 1/1
    variant_info_table[1].genotype_likelihoods.increment_by_index(4, 0.3L); // 4 -> 1/2
    variant_info_table[1].genotype_likelihoods.increment_by_index(5, 0.2L); // 5 -> 2/2
    std::vector<uint32_t> selected_genotype_indices = variant_info_table[1].genotype_likelihoods.select_genotypes();
	variant_info_table[1].update_active_alleles(2, selected_genotype_indices);
    read_set->setEntryAlleles(variant_info_table[1].position, variant_info_table[1].active_alleles);
    read_set->resetTags();

    /* Second round of phasing */
    {
        PhasingDPTable table = PhasingDPTable(read_set, &variant_info_table, false);
        std::vector<uint32_t> expected_cluster_ids = {0, 1, 1, 3, 3, 1, 3, 1, 3, 1, 10, 11, 11};
        std::vector<bool> expected_cluster_status = {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        std::vector<uint32_t> expected_constrained_cluster_ids = {-1u, 3, 3, 1, 1, 3, 1, 3, 1, 3, 11, 10, 10};
        std::vector<bool> expected_constrained_cluster_status = {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        assert(expected_cluster_ids.size() == read_set->size());
        for (uint32_t i = 0; i < read_set->size(); i++) {
            //std::cout << read_set->get(i)->getName() << "\t\t" << read_set->get(i)->isClustered() << "\t\t" << read_set->get(i)->getClusterID() <<"\t\t" << read_set->get(i)->getConstrainedClusterID() << std::endl;
            assert(expected_cluster_status[i] == read_set->get(i)->isClustered());
            assert(expected_cluster_ids[i] == read_set->get(i)->getClusterID());
            assert(expected_constrained_cluster_status[i] == read_set->get(i)->hasConstrainedCluster());
            if (expected_constrained_cluster_status[i]) {
                assert(expected_constrained_cluster_ids[i] == read_set->get(i)->getConstrainedClusterID());
            }
        }
    }
    
    /* Updating GenotypeLikelihoods for variant at 400 so that only allele 0 and 1 are chosen*/
    variant_info_table[3].genotype_likelihoods.increment_by_index(0, 0.5L); // 0 -> 0/0
    variant_info_table[3].genotype_likelihoods.increment_by_index(1, 0.3L); // 1 -> 0/1
    variant_info_table[3].genotype_likelihoods.increment_by_index(2, 0.2L); // 2 -> 1/1
    selected_genotype_indices = variant_info_table[3].genotype_likelihoods.select_genotypes();
	variant_info_table[3].update_active_alleles(2, selected_genotype_indices);
    read_set->setEntryAlleles(variant_info_table[3].position, variant_info_table[3].active_alleles);
    read_set->resetTags();

    /* Third round of phasing */
    {
        PhasingDPTable table = PhasingDPTable(read_set, &variant_info_table, false);
        std::vector<uint32_t> expected_cluster_ids = {0, 1, 1, 3, 3, 1, 3, 1, 3, 1, 10, 11, 11};
        std::vector<bool> expected_cluster_status = {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        std::vector<uint32_t> expected_constrained_cluster_ids = {-1u, 3, 3, 1, 1, 3, 1, 3, 1, 3, 11, 10, 10};
        std::vector<bool> expected_constrained_cluster_status = {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        assert(expected_cluster_ids.size() == read_set->size());
        for (uint32_t i = 0; i < read_set->size(); i++) {
            //std::cout << read_set->get(i)->getName() << "\t\t" << read_set->get(i)->isClustered() << "\t\t" << read_set->get(i)->getClusterID() <<"\t\t" << read_set->get(i)->getConstrainedClusterID() << std::endl;
            assert(expected_cluster_status[i] == read_set->get(i)->isClustered());
            assert(expected_cluster_ids[i] == read_set->get(i)->getClusterID());
            assert(expected_constrained_cluster_status[i] == read_set->get(i)->hasConstrainedCluster());
            if (expected_constrained_cluster_status[i]) {
                assert(expected_constrained_cluster_ids[i] == read_set->get(i)->getConstrainedClusterID());
            }
        }
        
    }

    delete read_set;

    assert_msg(true, "PhasingDPTable", "Multiposition DP table.");
}

void test_phasingdptable() {
    //test_singleposition_allreadsselected();
    //test_singleposition_somereadsselected();
    //test_singleposition_sv();
    test_multiposition();
}