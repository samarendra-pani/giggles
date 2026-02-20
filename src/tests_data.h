#ifndef TESTS_DATA_H
#define TESTS_DATA_H

#include <cassert>

#include "variantinfo.h"
#include "readset.h"

// Helper to print success messages
void assert_msg(bool condition, const std::string& prefix, const std::string& message) {
    if (!condition) {
        std::cerr << "[" << prefix << "] FAILED: " << message << std::endl;
        std::exit(1);
    }
    std::cout << "[" << prefix << "] PASSED: " << message << std::endl;
}

std::vector<variant_information_t> mock_variant_info_table_1() {

    std::vector<variant_information_t> variant_info_table;
    std::vector<uint32_t> position = {100, 200};
    uint32_t ploidy = 2;
    std::vector<uint32_t> n_alleles = {2, 2};
    std::vector<std::vector<int>> allele_references = { {1, -1, 1, 0}, 
                                                        {0, 1, 0, 1}};
    std::vector<bool> is_sv_position = {false, true};

    variant_info_table.push_back(variant_information_t(position[0], ploidy, n_alleles[0], allele_references[0], is_sv_position[0]));
    variant_info_table.push_back(variant_information_t(position[1], ploidy, n_alleles[1], allele_references[1], is_sv_position[1]));

    /**
     * Variant Information Table:
     * ........................| variant 1 | variant 2 
     * ------------------------|-----------|-----------
     * Position                | 100       | 200
     * Number of Alleles       | 2         | 2
     * Is a SV?                | No        | Yes
     * ------------------------------------------------
     * Haplotype 1 Allele      | 1         | 0
     * Haplotype 2 Allele      | -1        | 1
     * Haplotype 3 Allele      | 1         | 0
     * Haplotype 4 Allele      | 0         | 1
     */

    return variant_info_table;
}

ReadSet* mock_readset_1() {
    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setID(10);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setID(11);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setID(12);
    Read* read4 = new Read("read4", 60, 0); read_set->add(read4); read4->setID(13);
    Read* read5 = new Read("read5", 60, 0); read_set->add(read5); read5->setID(14);
    Read* read6 = new Read("read6", 60, 0); read_set->add(read6); read6->setID(15);
    Read* read7 = new Read("read7", 60, 0); read_set->add(read7); read7->setID(16);
    Read* read8 = new Read("read8", 60, 0); read_set->add(read8); read8->setID(17);
    Read* read9 = new Read("read9", 60, 0); read_set->add(read9); read9->setID(18);
    Read* read10 = new Read("read10", 60, 0); read_set->add(read10); read10->setID(19);
    Read* read11 = new Read("read11", 60, 0); read_set->add(read11); read11->setID(20);
    Read* read12 = new Read("read12", 60, 0); read_set->add(read12); read12->setID(21);
    Read* read13 = new Read("read13", 60, 0); read_set->add(read13); read13->setID(22);
    Read* read14 = new Read("read14", 60, 0); read_set->add(read14); read14->setID(23);
    Read* read15 = new Read("read15", 60, 0); read_set->add(read15); read15->setID(24);
    Read* read16 = new Read("read16", 60, 0); read_set->add(read16); read16->setID(25);
    Read* read17 = new Read("read17", 60, 0); read_set->add(read17); read17->setID(26);
    Read* read18 = new Read("read18", 60, 0); read_set->add(read18); read18->setID(27);
    Read* read19 = new Read("read19", 60, 0); read_set->add(read19); read19->setID(28);
    Read* read20 = new Read("read20", 60, 0); read_set->add(read20); read20->setID(29);
    
    
    /**
     * Adding variants to the reads
     */
    {
        std::vector<uint32_t> scores = std::vector<uint32_t>{10, 90};
        read1->addVariant(100, scores);
        read2->addVariant(100, scores);
        read3->addVariant(100, scores); read3->addVariant(200, scores);
        read4->addVariant(100, scores);
        read5->addVariant(100, scores); read5->addVariant(200, scores);
        read6->addVariant(100, scores);
        read7->addVariant(100, scores); read7->addVariant(200, scores);
        read8->addVariant(100, scores); read8->addVariant(200, scores);
        read9->addVariant(100, scores);
        read10->addVariant(100, scores); read10->addVariant(200, scores);
        read11->addVariant(100, scores); read11->addVariant(200, scores);
        read12->addVariant(100, scores); read12->addVariant(200, scores);
        read13->addVariant(100, scores); read13->addVariant(200, scores);
        read14->addVariant(100, scores);
        read15->addVariant(100, scores); read15->addVariant(200, scores);
        read16->addVariant(100, scores); read16->addVariant(200, scores);
        read17->addVariant(100, scores); read17->addVariant(200, scores);
        read18->addVariant(100, scores); read18->addVariant(200, scores);
        read19->addVariant(200, scores);
        read20->addVariant(200, scores);
    }
    /**
     * Setting up clusters
     */
    {
        read1->setClusterID(10); read1->setClusterStatus(true);
        read2->setClusterID(10); read2->setClusterStatus(true);

        read3->setClusterID(8); read3->setClusterStatus(true); read3->setConstrainedClusterID(1); 
        read7->setClusterID(8); read7->setClusterStatus(true); read7->setConstrainedClusterID(1);
        
        read4->setClusterID(7); read4->setClusterStatus(true); read4->setConstrainedClusterID(2);
        read5->setClusterID(7); read5->setClusterStatus(true); read5->setConstrainedClusterID(2);
        read14->setClusterID(7); read14->setClusterStatus(true); read14->setConstrainedClusterID(2);
        read15->setClusterID(2); read15->setClusterStatus(true); read15->setConstrainedClusterID(7);

        read8->setClusterID(17); read8->setClusterStatus(true); read8->setConstrainedClusterID(29);
        read19->setClusterID(17); read19->setClusterStatus(true); read19->setConstrainedClusterID(29);
        read20->setClusterID(29); read20->setClusterStatus(true); read20->setConstrainedClusterID(17);

        read9->setClusterID(3); read9->setClusterStatus(true);
        read10->setClusterID(3); read10->setClusterStatus(true);
        read13->setClusterID(3); read13->setClusterStatus(true);

        read11->setClusterID(9); read11->setClusterStatus(true); read11->setConstrainedClusterID(4);
        read12->setClusterID(9); read12->setClusterStatus(true); read12->setConstrainedClusterID(4);
        read17->setClusterID(9); read17->setClusterStatus(true); read17->setConstrainedClusterID(4);
        read16->setClusterID(4); read16->setClusterStatus(true); read16->setConstrainedClusterID(9);
        read18->setClusterID(4); read18->setClusterStatus(true); read18->setConstrainedClusterID(9);
        
    }

    /**
     * Reads Summary:
     * ID  | Name    | Variants       | ClusterID | Constraint
     * ----|---------|----------------|-----------|------------
     * 10  | read1   | 100            | 10        | -
     * 11  | read2   | 100            | 10        | -
     * 12  | read3   | 100, 200       | 8         | 1 
     * 13  | read4   | 100            | 7         | 2 
     * 14  | read5   | 100, 200       | 7         | 2 
     * 15  | read6   | 100            | -         | - 
     * 16  | read7   | 100, 200       | 8         | 1 
     * 17  | read8   | 100, 200       | 17        | 29 
     * 18  | read9   | 100            | 3         | - 
     * 19  | read10  | 100, 200       | 3         | - 
     * 20  | read11  | 100, 200       | 9         | 4 
     * 21  | read12  | 100, 200       | 9         | 4 
     * 22  | read13  | 100, 200       | 3         | - 
     * 23  | read14  | 100            | 7         | 2 
     * 24  | read15  | 100, 200       | 2         | 7 
     * 25  | read16  | 100, 200       | 4         | 9 
     * 26  | read17  | 100, 200       | 9         | 4 
     * 27  | read18  | 100, 200       | 4         | 9 
     * 28  | read19  | 200            | 17        | 29 
     * 29  | read20  | 200            | 29        | 17 
     */

    return read_set;
}

std::vector<variant_information_t> mock_variant_info_table_2() {

    std::vector<variant_information_t> variant_info_table;
    std::vector<uint32_t> position = {100, 200, 300, 400, 500};
    uint32_t ploidy = 2;
    std::vector<uint32_t> n_alleles = {2, 3, 2, 4, 2};
    std::vector<std::vector<int>> allele_references = { {1, 0, 1, 0}, 
                                                        {0, 1, 0, 2},
                                                        {1, 1, 0, 0},
                                                        {0, 1, 3, 2},
                                                        {0, 0, 1, 0}};
    std::vector<bool> is_sv_position = {false, true, false, true, false};

    variant_info_table.push_back(variant_information_t(position[0], ploidy, n_alleles[0], allele_references[0], is_sv_position[0]));
    variant_info_table.push_back(variant_information_t(position[1], ploidy, n_alleles[1], allele_references[1], is_sv_position[1]));
    variant_info_table.push_back(variant_information_t(position[2], ploidy, n_alleles[2], allele_references[2], is_sv_position[2]));
    variant_info_table.push_back(variant_information_t(position[3], ploidy, n_alleles[3], allele_references[3], is_sv_position[3]));
    variant_info_table.push_back(variant_information_t(position[4], ploidy, n_alleles[4], allele_references[4], is_sv_position[4]));

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
    return variant_info_table;
}

ReadSet* mock_readset_2() {
    
    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3);
    Read* read4 = new Read("read4", 60, 0); read_set->add(read4);

    /**
     * Adding variants to the reads
     */
    std::vector<uint32_t> scores_1 = std::vector<uint32_t>{10, 90};
    std::vector<uint32_t> scores_2 = std::vector<uint32_t>{20, 30, 50};
    std::vector<uint32_t> scores_3 = std::vector<uint32_t>{85, 15};
    std::vector<uint32_t> scores_4 = std::vector<uint32_t>{5, 25, 35, 35};
    std::vector<uint32_t> scores_5 = std::vector<uint32_t>{40, 60};
    
    read1->addVariant(100, scores_1); read1->addVariant(200, scores_2); read1->addVariant(300, scores_3); read1->addVariant(400, scores_4);
    read2->addVariant(100, scores_1); read1->addVariant(200, scores_2); read1->addVariant(300, scores_3); read1->addVariant(400, scores_4); read1->addVariant(500, scores_5);
    read3->addVariant(200, scores_2); read3->addVariant(300, scores_3); read3->addVariant(400, scores_4);
    read4->addVariant(300, scores_3); read4->addVariant(400, scores_4); read4->addVariant(500, scores_5);

    read_set->sort();
    read_set->reassignReadIds();

    assert(read1->getID() == 0);
    assert(read2->getID() == 1);
    assert(read3->getID() == 2);
    assert(read4->getID() == 3);

    /**
     * Reads Summary:
     * ID | Name   | Variants                 
     * ---|--------|--------------------------
     * 0  | read1  | 100, 200, 300, 400       
     * 1  | read2  | 100, 200, 300, 400, 500
     * 2  | read3  | 200, 300, 400
     * 3  | read4  | 300, 400, 500
     */

    return read_set;
}

#endif // TESTS_DATA_H