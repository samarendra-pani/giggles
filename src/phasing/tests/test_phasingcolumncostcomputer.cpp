#include "test_phasingcolumncostcomputer.h"

void test_unphasable_position() {
    /**
     * Setting up variant info table
     */
    std::vector<variant_information_t> variant_info_table;
    std::vector<uint32_t> position = {100};
    uint32_t ploidy = 2;
    std::vector<uint32_t> n_alleles = {2};
    std::vector<std::vector<int>> allele_references = {{1, 0, 1, 0}};
    std::vector<bool> is_sv_position = {false};
    variant_info_table.push_back(variant_information_t(position[0], ploidy, n_alleles[0], allele_references[0], is_sv_position[0]));
    variant_info_table[0].phasable = false;
    assert(variant_info_table[0].genotype_likelihoods.size() == 3);

    /**
     * Setting up reads
     */
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
    read2->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1
    read3->addVariant(100, std::vector<uint32_t>{10, 10}); // EQUAL_SCORES
    read4->addVariant(100, std::vector<uint32_t>{10, 10}); // EQUAL_SCORES
    read5->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1
    read6->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1
    read7->addVariant(100, std::vector<uint32_t>{20, 10}); // ALLELE2
    read8->addVariant(100, std::vector<uint32_t>{10, 20}); // ALLELE1
    read9->addVariant(100, std::vector<uint32_t>{5, 10}); // ALLELE1
    read10->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1

    PhasingColumnIterator* iterator = new PhasingColumnIterator(*read_set, &variant_info_table, false);
    std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next();
    PhasingColumnCostComputer* cost_computer = new PhasingColumnCostComputer(*column, variant_info_table.at(0));
    uint32_t cost;
    PhasingColumnCostComputer::phased_variant_t alleles;

    /** Checking all possible partitions */
    for (uint32_t i = 0; i < (1 << column->size()); i++) {
        cost_computer->set_partitioning(i);
        cost = cost_computer->get_cost();
        alleles = cost_computer->get_alleles();
        assert(cost == 0);
        assert(alleles.allele0 == Entry::EQUAL_SCORES);
        assert(alleles.allele1 == Entry::EQUAL_SCORES);
    }
    /** Checking for update function */
    for (uint32_t i = 0; i < column->size(); i++) {
        cost_computer->set_partitioning(0);
        cost_computer->update_partitioning(i);
        cost = cost_computer->get_cost();
        alleles = cost_computer->get_alleles();
        assert(cost == 0);
        assert(alleles.allele0 == Entry::EQUAL_SCORES);
        assert(alleles.allele1 == Entry::EQUAL_SCORES);
    }
    assert_msg(true, "PhasingColumnCostComputer", "Cost and alleles of unphasable position.");
}

void test_homozygous_position() {
    /**
     * Setting up variant info table
     */
    std::vector<variant_information_t> variant_info_table;
    std::vector<uint32_t> position = {100};
    uint32_t ploidy = 2;
    std::vector<uint32_t> n_alleles = {4};
    std::vector<std::vector<int>> allele_references = {{1, 2, 3, 0}};
    std::vector<bool> is_sv_position = {false};
    variant_info_table.push_back(variant_information_t(position[0], ploidy, n_alleles[0], allele_references[0], is_sv_position[0]));
    assert(variant_info_table[0].genotype_likelihoods.size() == 10);
    variant_info_table[0].genotype_likelihoods.set_by_index(0, 0.0L);   // 0/0
    variant_info_table[0].genotype_likelihoods.set_by_index(1, 0.0L);   // 0/1
    variant_info_table[0].genotype_likelihoods.set_by_index(2, 0.0L);   // 1/1
    variant_info_table[0].genotype_likelihoods.set_by_index(3, 0.02L);  // 0/2
    variant_info_table[0].genotype_likelihoods.set_by_index(4, 0.01L);  // 1/2
    variant_info_table[0].genotype_likelihoods.set_by_index(5, 0.95L);  // 2/2 (Active)
    variant_info_table[0].genotype_likelihoods.set_by_index(6, 0.0L);   // 0/3
    variant_info_table[0].genotype_likelihoods.set_by_index(7, 0.0L);   // 1/3
    variant_info_table[0].genotype_likelihoods.set_by_index(8, 0.02L);  // 2/3
    variant_info_table[0].genotype_likelihoods.set_by_index(9, 0.0L);   // 3/3
    std::vector<uint32_t> selected_genotype_indices = variant_info_table[0].genotype_likelihoods.select_genotypes();
    assert(selected_genotype_indices.size() == 1);
    assert(selected_genotype_indices[0] == 5);
    variant_info_table[0].update_active_alleles(2, selected_genotype_indices);
    assert(variant_info_table[0].count_active_alleles() == 1);

    /** Homozygous position determined by genotype likelihoods */
    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(true);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(true);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(true);
    Read* read4 = new Read("read4", 60, 0); read_set->add(read4); read4->setSelected(true);
    Read* read5 = new Read("read5", 60, 0); read_set->add(read5); read5->setSelected(true);

    read1->addVariant(100, std::vector<uint32_t>{10, 90, 50, 20});
    read2->addVariant(100, std::vector<uint32_t>{20, 40, 20, 10});
    read3->addVariant(100, std::vector<uint32_t>{10, 10, 30, 20});
    read4->addVariant(100, std::vector<uint32_t>{10, 30, 30, 10});
    read5->addVariant(100, std::vector<uint32_t>{20, 20, 30, 10});

    PhasingColumnIterator* iterator = new PhasingColumnIterator(*read_set, &variant_info_table, false);
    std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next();
    PhasingColumnCostComputer* cost_computer = new PhasingColumnCostComputer(*column, variant_info_table.at(0));
    uint32_t cost;
    PhasingColumnCostComputer::phased_variant_t alleles;

    /** Checking all possible partitions */
    for (uint32_t i = 0; i < (1 << column->size()); i++) {
        cost_computer->set_partitioning(i);
        cost = cost_computer->get_cost();
        alleles = cost_computer->get_alleles();
        assert(cost == 0);
        assert(alleles.allele0 == Entry::ALLELE1);
        assert(alleles.allele1 == Entry::ALLELE1);
    }
    for (uint32_t i = 0; i < column->size(); i++) {
        cost_computer->set_partitioning(0);
        cost_computer->update_partitioning(i);
        cost = cost_computer->get_cost();
        alleles = cost_computer->get_alleles();
        assert(cost == 0);
        assert(alleles.allele0 == Entry::ALLELE1);
        assert(alleles.allele1 == Entry::ALLELE1);
    }

    /** Homozygous position determined by pile-up */
    delete read_set;
    delete iterator;
    delete cost_computer;

    position = {100};
    ploidy = 2;
    n_alleles = {2};
    allele_references = {{1, 0, 1, 0}};
    is_sv_position = {false};
    variant_info_table[0] = variant_information_t(position[0], ploidy, n_alleles[0], allele_references[0], is_sv_position[0]);
    
    read_set = new ReadSet();
    read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(true);
    read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(true);
    read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(true);
    read4 = new Read("read4", 60, 0); read_set->add(read4); read4->setSelected(true);
    read5 = new Read("read5", 60, 0); read_set->add(read5); read5->setSelected(true);

    /** All reads point to ALLELE2 */
    read1->addVariant(100, std::vector<uint32_t>{10, 5});
    read2->addVariant(100, std::vector<uint32_t>{20, 2});
    read3->addVariant(100, std::vector<uint32_t>{10, 1});
    read4->addVariant(100, std::vector<uint32_t>{10, 5});
    read5->addVariant(100, std::vector<uint32_t>{20, 5});

    iterator = new PhasingColumnIterator(*read_set, &variant_info_table, false);
    column = iterator->get_next();
    cost_computer = new PhasingColumnCostComputer(*column, variant_info_table.at(0));
    
    /** Checking all possible partitions */
    for (uint32_t i = 0; i < (1 << column->size()); i++) {
        cost_computer->set_partitioning(i);
        cost = cost_computer->get_cost();
        alleles = cost_computer->get_alleles();
        assert(cost == 0);
        assert(alleles.allele0 == Entry::ALLELE2);
        assert(alleles.allele1 == Entry::ALLELE2);
    }
    for (uint32_t i = 0; i < column->size(); i++) {
        cost_computer->set_partitioning(0);
        cost_computer->update_partitioning(i);
        cost = cost_computer->get_cost();
        alleles = cost_computer->get_alleles();
        assert(cost == 0);
        assert(alleles.allele0 == Entry::ALLELE2);
        assert(alleles.allele1 == Entry::ALLELE2);
    }
    assert_msg(true, "PhasingColumnCostComputer", "Cost and alleles of homozygous position.");
}

void test_setting_partition_no_gl() {
    /**
     * Setting up variant info table
     */
    std::vector<variant_information_t> variant_info_table;
    std::vector<uint32_t> position = {100};
    uint32_t ploidy = 2;
    std::vector<uint32_t> n_alleles = {2};
    std::vector<std::vector<int>> allele_references = {{1, 0, 1, 0}};
    std::vector<bool> is_sv_position = {false};
    variant_info_table.push_back(variant_information_t(position[0], ploidy, n_alleles[0], allele_references[0], is_sv_position[0]));
    variant_info_table[0].phasable = true;
    assert(variant_info_table[0].genotype_likelihoods.size() == 3);

    /**
     * Setting up reads
     */
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
    read2->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1
    read3->addVariant(100, std::vector<uint32_t>{10, 10}); // EQUAL_SCORES
    read4->addVariant(100, std::vector<uint32_t>{10, 10}); // EQUAL_SCORES
    read5->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1
    read6->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1
    read7->addVariant(100, std::vector<uint32_t>{20, 10}); // ALLELE2
    read8->addVariant(100, std::vector<uint32_t>{10, 20}); // ALLELE1
    read9->addVariant(100, std::vector<uint32_t>{5, 10}); // ALLELE1
    read10->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1

    PhasingColumnIterator* iterator = new PhasingColumnIterator(*read_set, &variant_info_table, false);
    std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next();
    PhasingColumnCostComputer* cost_computer = new PhasingColumnCostComputer(*column, variant_info_table.at(0));
    uint32_t cost;
    PhasingColumnCostComputer::phased_variant_t alleles;


    /**
     * Partition 0 means all entries are in bipartition 0
     * Best cost is to flip read 7 to ALLELE1.
     */
    cost_computer->set_partitioning(0);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 30);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::EQUAL_SCORES);
    
    /**
     * In the above bipartition assignment of all reads being in the same bipartition,
     * the results should be unchanged if we send the BLANK and EQUAL_SCORES entry to the other bipartition.
     */
    cost_computer->set_partitioning(12);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 30);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::EQUAL_SCORES);
    
    
    /**
     * Partition 1 means all entries except read 1 are in bipartition 0 and read 1 is in bipartition 1.
     * The best cost is bipartition 1 in ALLELE1 and Bipartition 2 in ALLELE1.
     */
    cost_computer->set_partitioning(1);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 30);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::ALLELE1);
    
    /**
     * In the above bipartition assignment of all reads being in the same bipartition,
     * the results should be unchanged if we send the BLANK and EQUAL_SCORES entry to the other bipartition.
     */
    cost_computer->set_partitioning(13);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 30);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::ALLELE1);

     /**
     * Partition 64 means read 7 is in different bipartition.
     * Best cost is 0 with biparition 1 having ALLELE1 and bipartition 2 with ALLELE2
     */
    cost_computer->set_partitioning(64);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 0);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::ALLELE2);

    /**
     * In the above bipartition assignment of all reads being in the same bipartition,
     * the results should be unchanged if we send the BLANK and EQUAL_SCORES entry to the other bipartition.
     */
    cost_computer->set_partitioning(76);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 0);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::ALLELE2);

    assert_msg(true, "PhasingColumnCostComputer", "Setting partition with no prior genotype likelihoods.");
}

void test_setting_partition_with_gl() {
    /**
     * Setting up variant info table
     */
    std::vector<variant_information_t> variant_info_table;
    std::vector<uint32_t> position = {100};
    uint32_t ploidy = 2;
    std::vector<uint32_t> n_alleles = {2};
    std::vector<std::vector<int>> allele_references = {{1, 0, 1, 0}};
    std::vector<bool> is_sv_position = {false};
    variant_info_table.push_back(variant_information_t(position[0], ploidy, n_alleles[0], allele_references[0], is_sv_position[0]));
    assert(variant_info_table[0].genotype_likelihoods.size() == 3);
    /** Setting genotype likelihoods such that only HET is allowed for phasing. */
    variant_info_table[0].genotype_likelihoods.set_by_index(0, 0.05L);
    variant_info_table[0].genotype_likelihoods.set_by_index(1, 0.94L);
    variant_info_table[0].genotype_likelihoods.set_by_index(2, 0.01L);
    std::vector<uint32_t> selected_genotype_indices = variant_info_table[0].genotype_likelihoods.select_genotypes();
    assert(selected_genotype_indices.size() == 1);
    assert(selected_genotype_indices[0] == 1);
    variant_info_table[0].update_active_alleles(2, selected_genotype_indices);
    assert(variant_info_table[0].count_active_alleles() == 2);

    /**
     * Setting up reads
     */
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
    read2->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1
    read3->addVariant(100, std::vector<uint32_t>{10, 10}); // EQUAL_SCORES
    read4->addVariant(100, std::vector<uint32_t>{10, 10}); // EQUAL_SCORES
    read5->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1
    read6->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1
    read7->addVariant(100, std::vector<uint32_t>{20, 10}); // ALLELE2
    read8->addVariant(100, std::vector<uint32_t>{10, 20}); // ALLELE1
    read9->addVariant(100, std::vector<uint32_t>{5, 10}); // ALLELE1
    read10->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1

    PhasingColumnIterator* iterator = new PhasingColumnIterator(*read_set, &variant_info_table, false);
    std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next();
    PhasingColumnCostComputer* cost_computer = new PhasingColumnCostComputer(*column, variant_info_table.at(0));
    uint32_t cost;
    PhasingColumnCostComputer::phased_variant_t alleles;


    /**
     * Partition 0 means all entries are in bipartition 0
     * Best cost is to flip read 7 to ALLELE1.
     */
    cost_computer->set_partitioning(0);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 30);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::ALLELE2);
    
    /**
     * In the above bipartition assignment of all reads being in the same bipartition,
     * the results should be unchanged if we send the BLANK and EQUAL_SCORES entry to the other bipartition.
     */
    cost_computer->set_partitioning(12);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 30);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::ALLELE2);
    
    
    /**
     * Partition 1 means all entries except read 1 are in bipartition 0 and read 1 is in bipartition 1.
     * The best cost is bipartition 1 in ALLELE1 and Bipartition 2 in ALLELE2.
     */
    cost_computer->set_partitioning(1);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 60);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::ALLELE2);
    
    /**
     * In the above bipartition assignment of all reads being in the same bipartition,
     * the results should be unchanged if we send the BLANK and EQUAL_SCORES entry to the other bipartition.
     */
    cost_computer->set_partitioning(13);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 60);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::ALLELE2);

     /**
     * Partition 64 means read 7 is in different bipartition.
     * Best cost is 0 with biparition 1 having ALLELE1 and bipartition 2 with ALLELE2
     */
    cost_computer->set_partitioning(64);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 0);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::ALLELE2);

    /**
     * In the above bipartition assignment of all reads being in the same bipartition,
     * the results should be unchanged if we send the BLANK and EQUAL_SCORES entry to the other bipartition.
     */
    cost_computer->set_partitioning(76);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 0);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::ALLELE2);

    assert_msg(true, "PhasingColumnCostComputer", "Setting partition with prior genotype likelihoods.");
}

void test_update_partition_no_gl() {
    /**
     * Setting up variant info table
     */
    std::vector<variant_information_t> variant_info_table;
    std::vector<uint32_t> position = {100};
    uint32_t ploidy = 2;
    std::vector<uint32_t> n_alleles = {2};
    std::vector<std::vector<int>> allele_references = {{1, 0, 1, 0}};
    std::vector<bool> is_sv_position = {false};
    variant_info_table.push_back(variant_information_t(position[0], ploidy, n_alleles[0], allele_references[0], is_sv_position[0]));
    variant_info_table[0].phasable = true;
    assert(variant_info_table[0].genotype_likelihoods.size() == 3);

    /**
     * Setting up reads
     */
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
    read2->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1
    read3->addVariant(100, std::vector<uint32_t>{10, 10}); // EQUAL_SCORES
    read4->addVariant(100, std::vector<uint32_t>{10, 10}); // EQUAL_SCORES
    read5->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1
    read6->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1
    read7->addVariant(100, std::vector<uint32_t>{20, 10}); // ALLELE2
    read8->addVariant(100, std::vector<uint32_t>{10, 20}); // ALLELE1
    read9->addVariant(100, std::vector<uint32_t>{5, 10}); // ALLELE1
    read10->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1

    PhasingColumnIterator* iterator = new PhasingColumnIterator(*read_set, &variant_info_table, false);
    std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next();
    PhasingColumnCostComputer* cost_computer = new PhasingColumnCostComputer(*column, variant_info_table.at(0));
    uint32_t cost;
    PhasingColumnCostComputer::phased_variant_t alleles;


    cost_computer->set_partitioning(0);
    cost_computer->update_partitioning(0);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 30);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::ALLELE1);
    
    /**
     * Updating partition from 1 to 3 by flipping 1-index bit.
     */
    cost_computer->update_partitioning(1);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 30);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::ALLELE1);
    
    cost_computer->set_partitioning(0);
    alleles = cost_computer->get_alleles();
    cost_computer->update_partitioning(6);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 0);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::ALLELE2);

    cost_computer->update_partitioning(2);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 0);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::ALLELE2);

    cost_computer->update_partitioning(3);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 0);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::ALLELE2);

    cost_computer->update_partitioning(1);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 30);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::EQUAL_SCORES);

    cost_computer->update_partitioning(0);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 30);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::ALLELE1);

    assert_msg(true, "PhasingColumnCostComputer", "Updating partition with no prior genotype likelihoods.");
}

void test_update_partition_with_gl() {
    /**
     * Setting up variant info table
     */
    std::vector<variant_information_t> variant_info_table;
    std::vector<uint32_t> position = {100};
    uint32_t ploidy = 2;
    std::vector<uint32_t> n_alleles = {2};
    std::vector<std::vector<int>> allele_references = {{1, 0, 1, 0}};
    std::vector<bool> is_sv_position = {false};
    variant_info_table.push_back(variant_information_t(position[0], ploidy, n_alleles[0], allele_references[0], is_sv_position[0]));
    variant_info_table[0].phasable = true;
    assert(variant_info_table[0].genotype_likelihoods.size() == 3);
    /** Setting genotype likelihoods such that only HET is allowed for phasing. */
    variant_info_table[0].genotype_likelihoods.set_by_index(0, 0.12L);
    variant_info_table[0].genotype_likelihoods.set_by_index(1, 0.85L);
    variant_info_table[0].genotype_likelihoods.set_by_index(2, 0.03L);
    std::vector<uint32_t> selected_genotype_indices = variant_info_table[0].genotype_likelihoods.select_genotypes();
    assert(selected_genotype_indices.size() == 2);
    assert(selected_genotype_indices[0] == 1);
    assert(selected_genotype_indices[1] == 0);
    variant_info_table[0].update_active_alleles(2, selected_genotype_indices);
    assert(variant_info_table[0].count_active_alleles() == 2);

    /**
     * Setting up reads
     */
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
    read2->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1
    read3->addVariant(100, std::vector<uint32_t>{10, 10}); // EQUAL_SCORES
    read4->addVariant(100, std::vector<uint32_t>{10, 10}); // EQUAL_SCORES
    read5->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1
    read6->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1
    read7->addVariant(100, std::vector<uint32_t>{20, 10}); // ALLELE2
    read8->addVariant(100, std::vector<uint32_t>{10, 20}); // ALLELE1
    read9->addVariant(100, std::vector<uint32_t>{5, 10}); // ALLELE1
    read10->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1

    PhasingColumnIterator* iterator = new PhasingColumnIterator(*read_set, &variant_info_table, false);
    std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next();
    PhasingColumnCostComputer* cost_computer = new PhasingColumnCostComputer(*column, variant_info_table.at(0));
    uint32_t cost;
    PhasingColumnCostComputer::phased_variant_t alleles;

    cost_computer->set_partitioning(0);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 30);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::EQUAL_SCORES);
    
    cost_computer->update_partitioning(2);
    cost_computer->update_partitioning(3);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 30);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::EQUAL_SCORES);
    
    
    cost_computer->set_partitioning(0);
    cost_computer->update_partitioning(6);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 0);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::ALLELE2);

    cost_computer->update_partitioning(2);
    cost_computer->update_partitioning(3);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 0);
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::ALLELE2);

    assert_msg(true, "PhasingColumnCostComputer", "Updating partition with prior genotype likelihoods.");

}

void test_equal_scores() {
    /**
     * Setting up variant info table
     */
    std::vector<variant_information_t> variant_info_table;
    std::vector<uint32_t> position = {100};
    uint32_t ploidy = 2;
    std::vector<uint32_t> n_alleles = {2};
    std::vector<std::vector<int>> allele_references = {{1, 0, 1, 0}};
    std::vector<bool> is_sv_position = {false};
    variant_info_table.push_back(variant_information_t(position[0], ploidy, n_alleles[0], allele_references[0], is_sv_position[0]));
    assert(variant_info_table[0].genotype_likelihoods.size() == 3);

    /**
     * Setting up reads
     */
    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(true);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(true);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(true);
    Read* read4 = new Read("read4", 60, 0); read_set->add(read4); read4->setSelected(true);
    Read* read5 = new Read("read5", 60, 0); read_set->add(read5); read5->setSelected(true);

    read1->addVariant(100, std::vector<uint32_t>{10, 90}); // ALLELE1
    read2->addVariant(100, std::vector<uint32_t>{20, 10}); // ALLELE2
    read3->addVariant(100, std::vector<uint32_t>{10, 10}); // EQUAL_SCORES
    read4->addVariant(100, std::vector<uint32_t>{10, 10}); // EQUAL_SCORES
    read5->addVariant(100, std::vector<uint32_t>{20, 10}); // ALLELE2

    PhasingColumnIterator* iterator = new PhasingColumnIterator(*read_set, &variant_info_table, false);
    std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next();
    PhasingColumnCostComputer* cost_computer = new PhasingColumnCostComputer(*column, variant_info_table.at(0));
    uint32_t cost;
    PhasingColumnCostComputer::phased_variant_t alleles;

    cost_computer->set_partitioning(0);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 30);
    assert(alleles.allele0 == Entry::ALLELE2);
    assert(alleles.allele1 == Entry::EQUAL_SCORES);
    
    cost_computer->set_partitioning(2);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert(cost == 30);
    assert(alleles.allele0 == Entry::EQUAL_SCORES);
    assert(alleles.allele1 == Entry::ALLELE2);

    assert_msg(true, "PhasingColumnCostComputer", "Additional EQUAL_SCORE cases.");
}

void test_phasingcolumncostcomputer() {

    test_unphasable_position();
    test_homozygous_position();
    test_setting_partition_no_gl();
    test_setting_partition_with_gl();
    test_update_partition_no_gl();
    test_update_partition_with_gl();
    test_equal_scores();
}