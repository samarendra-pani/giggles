#ifndef TEST_PHASINGCOLUMNCOSTCOMPUTER_H
#define TEST_PHASINGCOLUMNCOSTCOMPUTER_H

#include "../phasingcolumncostcomputer.h"
#include "../phasingcolumniterator.h"
#include "../../tests_data.h"
#include <cassert>


void test_setting_partition() {
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
     * Based on the below genotype likelihoods, the PhredScore returned by
     * GenotypeLikelihoods::getPhredScore() are
     * 
     * 0 -> 0
     * 1 -> 5
     * 2 -> 8
     */
    variant_info_table[0].genotype_likelihoods.set_by_index(0, 0.7L);
    variant_info_table[0].genotype_likelihoods.set_by_index(1, 0.2L);
    variant_info_table[0].genotype_likelihoods.set_by_index(2, 0.1L);

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
    PhasingColumnCostComputer* cost_computer = new PhasingColumnCostComputer(*column, 0, &variant_info_table);
    uint32_t cost;
    PhasingColumnCostComputer::phased_variant_t alleles;


    /**
     * Partition 0 means all entries are in bipartition 0
     * The best cost is to flip read1 to ALLELE2 -> so Bipartition 0 has ALLELE2
     * Since no reads exist in Bipartition 1, Bipartition 1 can have ALLELE1 or ALLELE2.
     * 
     * Considering the Phred Scores, the best cost is if Bipartition 1 has ALLELE 1 since heterogygous genotype has lower score.
     * 
     * So best cost = 30 (flipping entry 1) + 5 (Phred Cost)
     * and best alleles are (ALLELE 2 | ALLELE 1).
     */
    cost_computer->set_partitioning(0);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert_msg(cost == 35, "PhasingColumnCostComputer", "Best cost for partition 0.");
    assert(alleles.allele0 == Entry::ALLELE2);
    assert(alleles.allele1 == Entry::ALLELE1);
    assert_msg(true, "PhasingColumnCostComputer", "Best alleles for partition 0.");
    /**
     * In the above bipartition assignment of all reads being in the same bipartition,
     * the results should be unchanged if we send the BLANK and EQUAL_SCORES entry to the other bipartition.
     */
    cost_computer->set_partitioning(12);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert_msg(cost == 35, "PhasingColumnCostComputer", "Best cost for partition 12 (Same as 0 but changed bipartition for EQUAL_SCORES and BLANK entries).");
    assert(alleles.allele0 == Entry::ALLELE2);
    assert(alleles.allele1 == Entry::ALLELE1);
    assert_msg(true, "PhasingColumnCostComputer", "Best alleles for partition 12 (Same as 0 but changed bipartition for EQUAL_SCORES and BLANK entries).");
    
    /**
     * Partition 1 means all entries except read 1 are in bipartition 0 and read 1 is in bipartition 1.
     * The best cost does not need any flipping and bipartition 0 has ALLELE2 and bipartition has ALLELE1.
     * 
     * So best cost = 5 (Phred Cost)
     * and best alleles are (ALLELE 2 | ALLELE 1).
     */
    cost_computer->set_partitioning(1);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert_msg(cost == 5, "PhasingColumnCostComputer", "Best cost for partition 1.");
    assert(alleles.allele0 == Entry::ALLELE2);
    assert(alleles.allele1 == Entry::ALLELE1);
    assert_msg(true, "PhasingColumnCostComputer", "Best alleles for partition 1.");
    /**
     * In the above bipartition assignment of all reads being in the same bipartition,
     * the results should be unchanged if we send the BLANK and EQUAL_SCORES entry to the other bipartition.
     */
    cost_computer->set_partitioning(13);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert_msg(cost == 5, "PhasingColumnCostComputer", "Best cost for partition 13 (Same as 1 but changed bipartition for EQUAL_SCORES and BLANK entries).");
    assert(alleles.allele0 == Entry::ALLELE2);
    assert(alleles.allele1 == Entry::ALLELE1);
    assert_msg(true, "PhasingColumnCostComputer", "Best alleles for partition 13 (Same as 1 but changed bipartition for EQUAL_SCORES and BLANK entries).");

     /**
     * Partition 18 means read 1, 3, 4 in bip0 and read 2, 5 in bip 1
     * The best cost does not need any flipping and bipartition 0 has ALLELE1 and bipartition has ALLELE2.
     * 
     * So best cost = 5 (Phred Cost)
     * and best alleles are (ALLELE 1 | ALLELE 2).
     */
    cost_computer->set_partitioning(18);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert_msg(cost == 5, "PhasingColumnCostComputer", "Best cost for partition 18.");
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::ALLELE2);
    assert_msg(true, "PhasingColumnCostComputer", "Best alleles for partition 18.");
    /**
     * In the above bipartition assignment of all reads being in the same bipartition,
     * the results should be unchanged if we send the BLANK and EQUAL_SCORES entry to the other bipartition.
     */
    cost_computer->set_partitioning(30);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert_msg(cost == 5, "PhasingColumnCostComputer", "Best cost for partition 30 (Same as 18 but changed bipartition for EQUAL_SCORES and BLANK entries).");
    assert(alleles.allele0 == Entry::ALLELE1);
    assert(alleles.allele1 == Entry::ALLELE2);
    assert_msg(true, "PhasingColumnCostComputer", "Best alleles for partition 30 (Same as 18 but changed bipartition for EQUAL_SCORES and BLANK entries).");
}

void test_update_partition() {
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
     * Based on the below genotype likelihoods, the PhredScore returned by
     * GenotypeLikelihoods::getPhredScore() are
     * 
     * 0 -> 0
     * 1 -> 5
     * 2 -> 8
     */
    variant_info_table[0].genotype_likelihoods.set_by_index(0, 0.7L);
    variant_info_table[0].genotype_likelihoods.set_by_index(1, 0.2L);
    variant_info_table[0].genotype_likelihoods.set_by_index(2, 0.1L);

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
    PhasingColumnCostComputer* cost_computer = new PhasingColumnCostComputer(*column, 0, &variant_info_table);
    uint32_t cost;
    PhasingColumnCostComputer::phased_variant_t alleles;

    /**
     * Partition 0 has best cost 35.
     * We flip the 0-index bit and make the partition 1
     */
    cost_computer->set_partitioning(0);
    cost_computer->update_partitioning(0);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert_msg(cost == 5, "PhasingColumnCostComputer", "Best cost for updating partition to 1.");
    assert(alleles.allele0 == Entry::ALLELE2);
    assert(alleles.allele1 == Entry::ALLELE1);
    assert_msg(true, "PhasingColumnCostComputer", "Best alleles for updating partition to 1.");

    /**
     * Updating partition from 1 to 3 by flipping 1-index bit.
     */
    cost_computer->update_partitioning(1);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert_msg(cost == 35, "PhasingColumnCostComputer", "Best cost for updating partition to 3.");
    assert(alleles.allele0 == Entry::ALLELE2);
    assert(alleles.allele1 == Entry::ALLELE1);
    assert_msg(true, "PhasingColumnCostComputer", "Best alleles for updating partition to 3.");
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
     * Based on the below genotype likelihoods, the PhredScore returned by
     * GenotypeLikelihoods::getPhredScore() are
     * 
     * 0 -> 1
     * 1 -> 0
     * 2 -> 1
     */
    variant_info_table[0].genotype_likelihoods.set_by_index(0, 0.3L);
    variant_info_table[0].genotype_likelihoods.set_by_index(1, 0.4L);
    variant_info_table[0].genotype_likelihoods.set_by_index(2, 0.3L);

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
    read3->addVariant(100, std::vector<uint32_t>{10, 20}); // ALLELE1
    read4->addVariant(100, std::vector<uint32_t>{10, 10}); // EQUAL_SCORES
    read5->addVariant(100, std::vector<uint32_t>{20, 10}); // ALLELE2

    PhasingColumnIterator* iterator = new PhasingColumnIterator(*read_set, &variant_info_table, false);
    std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next();
    PhasingColumnCostComputer* cost_computer = new PhasingColumnCostComputer(*column, 0, &variant_info_table);
    uint32_t cost;
    PhasingColumnCostComputer::phased_variant_t alleles;

    /**
     * Partition 0 has all entries in one partition.
     * 
     * But the best allele assignment can either be ALLELE1 | ALLELE2 or ALLELE2 | ALLELE1.
     * 
     * So best cost is 60 but best alleles can be multiple
     */
    cost_computer->set_partitioning(0);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert_msg(cost == 60, "PhasingColumnCostComputer", "Best cost for partition 0 with EQUAL SCORES.");
    assert(alleles.allele0 == Entry::EQUAL_SCORES);
    assert(alleles.allele1 == Entry::EQUAL_SCORES);
    assert_msg(true, "PhasingColumnCostComputer", "Best alleles for partition 0 with EQUAL SCORES.");
    /**
     * Should be the same for partition = 8
     */
    cost_computer->set_partitioning(8);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert_msg(cost == 60, "PhasingColumnCostComputer", "Best cost for partition 8 with EQUAL SCORES (same as partition 0).");
    assert(alleles.allele0 == Entry::EQUAL_SCORES);
    assert(alleles.allele1 == Entry::EQUAL_SCORES);
    assert_msg(true, "PhasingColumnCostComputer", "Best alleles for partition 8 with EQUAL SCORES (same as partition 0).");
}

void test_multiallelic_site() {
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
    /**
     * Based on the below genotype likelihoods, the PhredScore returned by
     * GenotypeLikelihoods::getPhredScore() are
     * 
     * 0 -> 7 (Active)
     * 1 -> 4
     * 2 -> 7
     * 3 -> 1
     * 4 -> 7
     * 5 -> 10
     * 6 -> 0 (Active)
     * 7 -> 4
     * 8 -> 4
     * 9 -> 10 (Active)
     */
    variant_info_table[0].genotype_likelihoods.set_by_index(0, 0.05L);  // 0/0 (Active)
    variant_info_table[0].genotype_likelihoods.set_by_index(1, 0.1L);   // 0/1
    variant_info_table[0].genotype_likelihoods.set_by_index(2, 0.05L);  // 1/1
    variant_info_table[0].genotype_likelihoods.set_by_index(3, 0.2L);   // 0/2
    variant_info_table[0].genotype_likelihoods.set_by_index(4, 0.05L);  // 1/2
    variant_info_table[0].genotype_likelihoods.set_by_index(5, 0.025L); // 2/2
    variant_info_table[0].genotype_likelihoods.set_by_index(6, 0.3L);   // 0/3 (Active)
    variant_info_table[0].genotype_likelihoods.set_by_index(7, 0.1L);   // 1/3
    variant_info_table[0].genotype_likelihoods.set_by_index(8, 0.1L);   // 2/3
    variant_info_table[0].genotype_likelihoods.set_by_index(9, 0.025L); // 3/3 (Active)
    variant_info_table[0].set_allele_inactive(1);
    variant_info_table[0].set_allele_inactive(2);

    /**
     * Setting up reads
     */
    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(true);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(true);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(true);
    Read* read4 = new Read("read4", 60, 0); read_set->add(read4); read4->setSelected(true);
    Read* read5 = new Read("read5", 60, 0); read_set->add(read5); read5->setSelected(true);

    read1->addVariant(100, std::vector<uint32_t>{10, 90, 50, 20}); // ALLELE1
    read2->addVariant(100, std::vector<uint32_t>{20, 40, 20, 10}); // ALLELE2
    read3->addVariant(100, std::vector<uint32_t>{10, 10, 30, 20}); // ALLELE1
    read4->addVariant(100, std::vector<uint32_t>{10, 30, 30, 10}); // EQUAL_SCORES
    read5->addVariant(100, std::vector<uint32_t>{20, 20, 30, 10}); // ALLELE2

    PhasingColumnIterator* iterator = new PhasingColumnIterator(*read_set, &variant_info_table, false);
    std::unique_ptr<std::vector<const Entry*> > column = iterator->get_next();
    PhasingColumnCostComputer* cost_computer = new PhasingColumnCostComputer(*column, 0, &variant_info_table);
    uint32_t cost;
    PhasingColumnCostComputer::phased_variant_t alleles;

    /**
     * Partition 0 means all entries are in bipartition 0
     * The best cost is to flip either read 1, 3 or read 2, 5.
     * 
     * Flip cost is 60 and Phred cost is 0. But (ALLELE 1 | ALLELE 2) (ALLELE2 | ALLELE1) is possible
     */
    cost_computer->set_partitioning(0);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert_msg(cost == 60, "PhasingColumnCostComputer", "Multiallelic: Best cost for partition 0.");
    assert(alleles.allele0 == Entry::EQUAL_SCORES);
    assert(alleles.allele1 == Entry::EQUAL_SCORES);
    assert_msg(true, "PhasingColumnCostComputer", "Multiallelic: Best alleles for partition 0.");
    /**
     * In the above bipartition assignment of all reads being in the same bipartition,
     * the results should be unchanged if we send the BLANK and EQUAL_SCORES entry to the other bipartition.
     */
    cost_computer->set_partitioning(8);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert_msg(cost == 60, "PhasingColumnCostComputer", "Multiallelic: Best cost for partition 8 (Same as 0 but changed bipartition for EQUAL_SCORES and BLANK entries).");
    assert(alleles.allele0 == Entry::EQUAL_SCORES);
    assert(alleles.allele1 == Entry::EQUAL_SCORES);
    assert_msg(true, "PhasingColumnCostComputer", "Multiallelic: Best alleles for partition 8 (Same as 0 but changed bipartition for EQUAL_SCORES and BLANK entries).");
    

    cost_computer->set_partitioning(0);
    cost_computer->update_partitioning(0); // partition is now 1
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert_msg(cost == 30, "PhasingColumnCostComputer", "Multiallelic: Best cost for partition 1.");
    assert(alleles.allele0 == Entry::ALLELE2);
    assert(alleles.allele1 == Entry::ALLELE1);
    assert_msg(true, "PhasingColumnCostComputer", "Multiallelic: Best alleles for partition 1.");
    
    cost_computer->set_partitioning(9);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert_msg(cost == 30, "PhasingColumnCostComputer", "Multiallelic: Best cost for partition 9 (Same as 1 but changed bipartition for EQUAL_SCORES and BLANK entries).");
    assert(alleles.allele0 == Entry::ALLELE2);
    assert(alleles.allele1 == Entry::ALLELE1);
    assert_msg(true, "PhasingColumnCostComputer", "Multiallelic: Best alleles for partition 9 (Same as 1 but changed bipartition for EQUAL_SCORES and BLANK entries).");

    
    
    cost_computer->set_partitioning(5);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert_msg(cost == 0, "PhasingColumnCostComputer", "Multiallelic: Best cost for partition 5.");
    assert(alleles.allele0 == Entry::ALLELE2);
    assert(alleles.allele1 == Entry::ALLELE1);
    assert_msg(true, "PhasingColumnCostComputer", "Multiallelic: Best alleles for partition 5.");

    cost_computer->set_partitioning(13);
    cost = cost_computer->get_cost();
    alleles = cost_computer->get_alleles();
    assert_msg(cost == 0, "PhasingColumnCostComputer", "Multiallelic: Best cost for partition 13 (Same as 5 but changed bipartition for EQUAL_SCORES and BLANK entries).");
    assert(alleles.allele0 == Entry::ALLELE2);
    assert(alleles.allele1 == Entry::ALLELE1);
    assert_msg(true, "PhasingColumnCostComputer", "Multiallelic: Best alleles for partition 13 (Same as 5 but changed bipartition for EQUAL_SCORES and BLANK entries).");
}

/**
 * TODO: More tests for cases where only haplotype has EQUAL SCORES
 */
void test_phasingcolumncostcomputer() {

    test_setting_partition();
    test_update_partition();
    test_equal_scores();
    test_multiallelic_site();
}

#endif // TEST_PHASINGCOLUMNCOSTCOMPUTER_H
