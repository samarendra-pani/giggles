#include "test_haplotypesampler.h"

void test_singlecolumn() {
    
    /**
     * Setting up variant info table
     */
    std::vector<variant_information_t> variant_info_table;
    std::vector<uint32_t> position = {100};
    uint32_t ploidy = 2;
    std::vector<uint32_t> n_alleles = {4};
    std::vector<std::vector<int>> allele_references = {{0, 1, 0, 1, 2, 0, 1, 3, 0, 1, 3, 2, -1, 0}};
    std::vector<bool> is_sv_position = {false};
    variant_info_table.push_back(variant_information_t(position[0], ploidy, n_alleles[0], allele_references[0], is_sv_position[0]));

    /**
     * Setting up reads
     */
    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(true);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(false);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(true);
    Read* read4 = new Read("read4", 60, 0); read_set->add(read4); read4->setSelected(false);
    Read* read5 = new Read("read5", 60, 0); read_set->add(read5); read5->setSelected(true);
    Read* read6 = new Read("read6", 60, 0); read_set->add(read6); read6->setSelected(false);
    Read* read7 = new Read("read7", 60, 0); read_set->add(read7); read7->setSelected(false);
    Read* read8 = new Read("read8", 60, 0); read_set->add(read8); read8->setSelected(true);
    Read* read9 = new Read("read9", 60, 0); read_set->add(read9); read9->setSelected(true);
    Read* read10 = new Read("read10", 60, 0); read_set->add(read10); read10->setSelected(true);

    read1->addVariant(100, std::vector<uint8_t>{90, 10, 0, 10});
    read2->addVariant(100, std::vector<uint8_t>{90, 10, 0, 10});
    read3->addVariant(100, std::vector<uint8_t>{90, 80, 0, 10});
    read4->addVariant(100, std::vector<uint8_t>{90, 80, 0, 10});
    read5->addVariant(100, std::vector<uint8_t>{90, 10, 0, 10});
    read6->addVariant(100, std::vector<uint8_t>{90, 10, 0, 10});
    read7->addVariant(100, std::vector<uint8_t>{90, 80, 0, 10});
    read8->addVariant(100, std::vector<uint8_t>{90, 80, 0, 10});
    read9->addVariant(100, std::vector<uint8_t>{95, 90, 0, 10});
    read10->addVariant(100, std::vector<uint8_t>{90, 10, 0, 10});
    
    HaplotypeSampler* sampler = new HaplotypeSampler(read_set, &variant_info_table, 2, 1.0, 1.0, nullptr, true, 10);
    std::vector<variant_information_t> updated_table = sampler->get_updated_variant_table(2);

    assert(updated_table.size() == 1);
    assert(updated_table[0].count_active_alleles() == 1);
    assert(updated_table[0].get_active_positions()[0] == 0);

    sampler = nullptr;
    sampler = new HaplotypeSampler(read_set, &variant_info_table, 6, 1.0, 1.0, nullptr, true, 10);
    updated_table = sampler->get_updated_variant_table(2);
    assert(updated_table.size() == 1);
    assert(updated_table[0].count_active_alleles() == 2);
    assert(updated_table[0].get_active_positions()[0] == 0);
    assert(updated_table[0].get_active_positions()[1] == 1);

    sampler = nullptr;
    sampler = new HaplotypeSampler(read_set, &variant_info_table, 8, 1.0, 1.0, nullptr, true, 10);
    updated_table = sampler->get_updated_variant_table(2);
    assert(updated_table.size() == 1);
    assert(updated_table[0].count_active_alleles() == 3);
    assert(updated_table[0].get_active_positions()[0] == 0);
    assert(updated_table[0].get_active_positions()[1] == 1);
    assert(updated_table[0].get_active_positions()[2] == 2);

    assert_msg(true, "HaplotypeSampler", "Single column haplotype sampling.");
}

void test_multiplecolumns() {
    /**
     * Setting up variant info table
     */
    std::vector<variant_information_t> variant_info_table;
    std::vector<uint32_t> position = {100, 200, 300};
    uint32_t ploidy = 2;
    std::vector<uint32_t> n_alleles = {4, 2, 3};
    std::vector<std::vector<int>> allele_references = {
        {0, 1, 0, 1, 2, 0, 1, 3, 0, 1, 3, 2,-1, 0},
        {0, 0, 0, 1, 0, 1, 0,-1, 0,-1, 1, 0, 1, 0},
        {0, 0, 1, 2, 0, 1, 2, 0, 0, 1, 1, 0, 2, 0}
    };
    std::vector<bool> is_sv_position = {false, false, false};
    variant_info_table.push_back(variant_information_t(position[0], ploidy, n_alleles[0], allele_references[0], is_sv_position[0]));
    variant_info_table.push_back(variant_information_t(position[1], ploidy, n_alleles[1], allele_references[1], is_sv_position[1]));
    variant_info_table.push_back(variant_information_t(position[2], ploidy, n_alleles[2], allele_references[2], is_sv_position[2]));

    /**
     * Setting up reads
     */
    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3);
    Read* read4 = new Read("read4", 60, 0); read_set->add(read4);
    Read* read5 = new Read("read5", 60, 0); read_set->add(read5);
    Read* read6 = new Read("read6", 60, 0); read_set->add(read6);
    Read* read7 = new Read("read7", 60, 0); read_set->add(read7);
    Read* read8 = new Read("read8", 60, 0); read_set->add(read8);
    Read* read9 = new Read("read9", 60, 0); read_set->add(read9);
    Read* read10 = new Read("read10", 60, 0); read_set->add(read10);
    {
        read1->addVariant(100, std::vector<uint8_t>{90, 10, 0, 10});
        read2->addVariant(100, std::vector<uint8_t>{10, 10, 90, 10});
        read3->addVariant(100, std::vector<uint8_t>{90, 10, 0, 10});
        read4->addVariant(100, std::vector<uint8_t>{10, 10, 90, 10});
        read5->addVariant(100, std::vector<uint8_t>{90, 10, 0, 10});
        read6->addVariant(100, std::vector<uint8_t>{10, 10, 90, 10});
        read7->addVariant(100, std::vector<uint8_t>{90, 10, 0, 10});
        read8->addVariant(100, std::vector<uint8_t>{10, 10, 90, 10});
        read9->addVariant(100, std::vector<uint8_t>{90, 10, 0, 10});
        read10->addVariant(100, std::vector<uint8_t>{90, 10, 0, 10});

        read1->addVariant(200, std::vector<uint8_t>{90, 10});
        read2->addVariant(200, std::vector<uint8_t>{10, 90});
        read3->addVariant(200, std::vector<uint8_t>{90, 10});
        read4->addVariant(200, std::vector<uint8_t>{10, 90});
        read5->addVariant(200, std::vector<uint8_t>{90, 10});
        read6->addVariant(200, std::vector<uint8_t>{10, 90});
        read7->addVariant(200, std::vector<uint8_t>{90, 10});
        read8->addVariant(200, std::vector<uint8_t>{10, 90});
        read9->addVariant(200, std::vector<uint8_t>{90, 10});
        read10->addVariant(200, std::vector<uint8_t>{90, 10});

        read1->addVariant(300, std::vector<uint8_t>{10, 90, 0});
        read2->addVariant(300, std::vector<uint8_t>{10, 0, 90});
        read3->addVariant(300, std::vector<uint8_t>{10, 90, 0});
        read4->addVariant(300, std::vector<uint8_t>{10, 0, 90});
        read5->addVariant(300, std::vector<uint8_t>{10, 90, 0});
        read6->addVariant(300, std::vector<uint8_t>{10, 0, 90});
        read7->addVariant(300, std::vector<uint8_t>{10, 90, 0});
        read8->addVariant(300, std::vector<uint8_t>{10, 0, 90});
        read9->addVariant(300, std::vector<uint8_t>{10, 90, 0});
        read10->addVariant(300, std::vector<uint8_t>{10, 90, 0});
    }

    std::vector<uint32_t>* best_scores = new std::vector<uint32_t>();
    HaplotypeSampler* sampler = new HaplotypeSampler(read_set, &variant_info_table, 4, 1.0, 25000.0F, best_scores, true, 10);
    
    std::vector<variant_information_t> updated_table = sampler->get_updated_variant_table(2);

    /*
    std::cout << "Best scores: ";
    for (auto b: *best_scores) {
        std::cout << b << " ";
    }
    std::cout << std::endl;
    */

    assert(updated_table.size() == 3);
    assert(updated_table[0].count_active_alleles() == 2);
    assert(updated_table[0].get_active_positions()[0] == 0);
    assert(updated_table[0].get_active_positions()[1] == 2);
    assert(updated_table[1].count_active_alleles() == 2);
    assert(updated_table[1].get_active_positions()[0] == 0);
    assert(updated_table[1].get_active_positions()[1] == 1);
    assert(updated_table[2].count_active_alleles() == 2);
    assert(updated_table[2].get_active_positions()[0] == 1);
    assert(updated_table[2].get_active_positions()[1] == 2);

    assert(best_scores->at(0) == 6);    // Emission for 0 (at 100) = 2, 0 (at 200) = 2, and 1 (at 300) = 2
    assert(best_scores->at(1) == 23);   // Emission for 2 (at 100) = 3, 1 (at 200) = 3, and 2 (at 300) = 3. And a haplotype change between 100 and 200 which is 14.
    assert(best_scores->at(2) == 37);   // Emission for 0 (at 100) = 12, 1 (at 200) = 13, and 1 (at 300) = 12
    assert(best_scores->at(3) == 52);   // Emission for 2 (at 100) = 13, 0 (at 200) = 12, and 2 (at 300) = 13. And a haplotype change which is 14.
    assert_msg(true, "HaplotypeSampler", "Multi-column haplotype sampling.");
}

void test_haplotypesampler() {
    test_singlecolumn();
    test_multiplecolumns();
}