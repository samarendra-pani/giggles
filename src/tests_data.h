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
    std::vector<uint32_t> position = {1, 100, 200};
    uint32_t ploidy = 2;
    std::vector<uint32_t> n_alleles = {2, 2, 2};
    std::vector<std::vector<int>> allele_references = { {1, 1, 1, 1},
                                                        {1, -1, 1, 0}, 
                                                        {0, 1, 0, 1}};
    std::vector<bool> is_sv_position = {false, false, true};

    variant_info_table.push_back(variant_information_t(position[0], ploidy, n_alleles[0], allele_references[0], is_sv_position[0]));
    variant_info_table.push_back(variant_information_t(position[1], ploidy, n_alleles[1], allele_references[1], is_sv_position[1]));
    variant_info_table.push_back(variant_information_t(position[2], ploidy, n_alleles[2], allele_references[2], is_sv_position[2]));

    /**
     * Variant Information Table:
     * ........................| variant 0 | variant 1 | variant 2 
     * ------------------------|-----------|----------------------
     * Position                | 1         | 100       | 200
     * Number of Alleles       | 2         | 2         | 2
     * Is a SV?                | No        | No        | Yes
     * -----------------------------------------------------------
     * Haplotype 1 Allele      | 1         | 1         | 0
     * Haplotype 2 Allele      | 1         | -1        | 1
     * Haplotype 3 Allele      | 1         | 1         | 0
     * Haplotype 4 Allele      | 1         | 0         | 1
     * 
     * NOTE: variant 0 is a dummy variant for the dummy reads
     */

    return variant_info_table;
}

ReadSet* mock_readset_1() {
    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); // dummy read
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); // dummy read
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3); // dummy read
    Read* read4 = new Read("read4", 60, 0); read_set->add(read4); // dummy read
    Read* read5 = new Read("read5", 60, 0); read_set->add(read5); // dummy read
    Read* read6 = new Read("read6", 60, 0); read_set->add(read6); // dummy read
    Read* read7 = new Read("read7", 60, 0); read_set->add(read7); // dummy read
    Read* read8 = new Read("read8", 60, 0); read_set->add(read8); // dummy read
    Read* read9 = new Read("read9", 60, 0); read_set->add(read9); // dummy read
    Read* read10 = new Read("read10", 60, 0); read_set->add(read10); // dummy read
    Read* read11 = new Read("read11", 60, 0); read_set->add(read11);
    Read* read12 = new Read("read12", 60, 0); read_set->add(read12);
    Read* read13 = new Read("read13", 60, 0); read_set->add(read13);
    Read* read14 = new Read("read14", 60, 0); read_set->add(read14);
    Read* read15 = new Read("read15", 60, 0); read_set->add(read15);
    Read* read16 = new Read("read16", 60, 0); read_set->add(read16);
    Read* read17 = new Read("read17", 60, 0); read_set->add(read17);
    Read* read18 = new Read("read18", 60, 0); read_set->add(read18);
    Read* read19 = new Read("read19", 60, 0); read_set->add(read19);
    Read* read20 = new Read("read20", 60, 0); read_set->add(read20);
    Read* read21 = new Read("read21", 60, 0); read_set->add(read21);
    Read* read22 = new Read("read22", 60, 0); read_set->add(read22);
    Read* read23 = new Read("read23", 60, 0); read_set->add(read23);
    Read* read24 = new Read("read24", 60, 0); read_set->add(read24);
    Read* read25 = new Read("read25", 60, 0); read_set->add(read25);
    Read* read26 = new Read("read26", 60, 0); read_set->add(read26);
    Read* read27 = new Read("read27", 60, 0); read_set->add(read27);
    Read* read28 = new Read("read28", 60, 0); read_set->add(read28);
    Read* read29 = new Read("read29", 60, 0); read_set->add(read29);
    Read* read30 = new Read("read30", 60, 0); read_set->add(read30);
    
    
    /**
     * Adding variants to the reads
     */
    {
        std::vector<uint32_t> scores = std::vector<uint32_t>{10, 90};
        read1->addVariant(1, scores);
        read2->addVariant(1, scores);
        read3->addVariant(1, scores);
        read4->addVariant(1, scores);
        read5->addVariant(1, scores);
        read6->addVariant(1, scores);
        read7->addVariant(1, scores);
        read8->addVariant(1, scores);
        read9->addVariant(1, scores);
        read10->addVariant(1, scores);
        read11->addVariant(100, scores);
        read12->addVariant(100, scores);
        read13->addVariant(100, scores); read13->addVariant(200, scores);
        read14->addVariant(100, scores);
        read15->addVariant(100, scores); read15->addVariant(200, scores);
        read16->addVariant(100, scores);
        read17->addVariant(100, scores); read17->addVariant(200, scores);
        read18->addVariant(100, scores); read18->addVariant(200, scores);
        read19->addVariant(100, scores);
        read20->addVariant(100, scores); read20->addVariant(200, scores);
        read21->addVariant(100, scores); read21->addVariant(200, scores);
        read22->addVariant(100, scores); read22->addVariant(200, scores);
        read23->addVariant(100, scores); read23->addVariant(200, scores);
        read24->addVariant(100, scores);
        read25->addVariant(100, scores); read25->addVariant(200, scores);
        read26->addVariant(100, scores); read26->addVariant(200, scores);
        read27->addVariant(100, scores); read27->addVariant(200, scores);
        read28->addVariant(100, scores); read28->addVariant(200, scores);
        read29->addVariant(200, scores);
        read30->addVariant(200, scores);
    }

    /**
     * Setting IDs
     */
    {   read1->setID(0);
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
        read14->setID(13);
        read15->setID(14);
        read16->setID(15);
        read17->setID(16);
        read18->setID(17);
        read19->setID(18);
        read20->setID(19);
        read21->setID(20);
        read22->setID(21);
        read23->setID(22);
        read24->setID(23);
        read25->setID(24);
        read26->setID(25);
        read27->setID(26);
        read28->setID(27);
        read29->setID(28);
        read30->setID(29);
    }


    /**
     * Setting up clusters
     */
    {
        read11->setClusterID(10); read11->setClusterStatus(true);
        read12->setClusterID(10); read12->setClusterStatus(true);

        read13->setClusterID(8); read13->setClusterStatus(true); read13->setConstrainedClusterID(1); 
        read17->setClusterID(8); read17->setClusterStatus(true); read17->setConstrainedClusterID(1);
        
        read14->setClusterID(7); read14->setClusterStatus(true); read14->setConstrainedClusterID(2);
        read15->setClusterID(7); read15->setClusterStatus(true); read15->setConstrainedClusterID(2);
        read24->setClusterID(7); read24->setClusterStatus(true); read24->setConstrainedClusterID(2);
        read25->setClusterID(2); read25->setClusterStatus(true); read25->setConstrainedClusterID(7);

        read18->setClusterID(17); read18->setClusterStatus(true); read18->setConstrainedClusterID(29);
        read29->setClusterID(17); read29->setClusterStatus(true); read29->setConstrainedClusterID(29);
        read30->setClusterID(29); read30->setClusterStatus(true); read30->setConstrainedClusterID(17);

        read19->setClusterID(3); read19->setClusterStatus(true);
        read20->setClusterID(3); read20->setClusterStatus(true);
        read23->setClusterID(3); read23->setClusterStatus(true);

        read21->setClusterID(9); read21->setClusterStatus(true); read21->setConstrainedClusterID(4);
        read22->setClusterID(9); read22->setClusterStatus(true); read22->setConstrainedClusterID(4);
        read27->setClusterID(9); read27->setClusterStatus(true); read27->setConstrainedClusterID(4);
        read26->setClusterID(4); read26->setClusterStatus(true); read26->setConstrainedClusterID(9);
        read28->setClusterID(4); read28->setClusterStatus(true); read28->setConstrainedClusterID(9);
        
    }

    /**
     * Reads Summary:
     * ID  | Name    | Variants       | ClusterID | Constraint
     * ----|---------|----------------|-----------|------------
     * 10  | read11  | 100            | 10        | -
     * 11  | read12  | 100            | 10        | -
     * 12  | read13  | 100, 200       | 8         | 1 
     * 13  | read14  | 100            | 7         | 2 
     * 14  | read15  | 100, 200       | 7         | 2 
     * 15  | read16  | 100            | -         | - 
     * 16  | read17  | 100, 200       | 8         | 1 
     * 17  | read18  | 100, 200       | 17        | 29 
     * 18  | read19  | 100            | 3         | - 
     * 19  | read20  | 100, 200       | 3         | - 
     * 20  | read21  | 100, 200       | 9         | 4 
     * 21  | read22  | 100, 200       | 9         | 4 
     * 22  | read23  | 100, 200       | 3         | - 
     * 23  | read24  | 100            | 7         | 2 
     * 24  | read25  | 100, 200       | 2         | 7 
     * 25  | read26  | 100, 200       | 4         | 9 
     * 26  | read27  | 100, 200       | 9         | 4 
     * 27  | read28  | 100, 200       | 4         | 9 
     * 28  | read29  | 200            | 17        | 29 
     * 29  | read30  | 200            | 29        | 17 
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