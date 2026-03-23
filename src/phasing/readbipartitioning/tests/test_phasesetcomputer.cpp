#include "test_phasesetcomputer.h"

void test_heterozygous_positions() {
    
    ReadSet* superreads = mock_superreads();
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_3();
    std::vector<uint32_t>* accessible_positions = new std::vector<uint32_t>();
	for (uint32_t i = 0; i < variant_info_table.size(); ++i) {
		if (variant_info_table[i].phasable) {
			accessible_positions->push_back(variant_info_table[i].position);
		}
	}
    assert(accessible_positions->size() == 13);
    std::unordered_set<uint32_t> heterozygous_positions;
    std::unordered_set<uint32_t> accessible_positions_set(accessible_positions->begin(), accessible_positions->end());

    assert (superreads->size() == 2); // two superreads represeting the two haplotypes
    Read* superread0 = superreads->get(0);
    Read* superread1 = superreads->get(1);

    assert (superread0->getVariantCount() == superread1->getVariantCount()); // both superreads should have same number of variants
    for (uint32_t i = 0; i < superread0->getVariantCount(); ++i) {
        assert (superread0->getPosition(i) == superread1->getPosition(i)); // both superreads should have variants at same positions
        // skip positions that are not accessible (that have more than 2 active alleles)
        if (accessible_positions_set.find(superread0->getPosition(i)) == accessible_positions_set.end()) { continue; }
        
        Entry::allele_t allele0 = superread0->getEntry(i)->get_allele_type();
        Entry::allele_t allele1 = superread1->getEntry(i)->get_allele_type();

        if ((allele0 == Entry::ALLELE1 && allele1 == Entry::ALLELE2) ||
            (allele0 == Entry::ALLELE2 && allele1 == Entry::ALLELE1)) {
            heterozygous_positions.insert(superread0->getPosition(i));
        }
    }
    assert(heterozygous_positions.size() == 7);
    std::vector<uint32_t> expected_het_positions = {100, 400, 600, 800, 900, 1100, 1500};
    std::unordered_set<uint32_t> tempSet(expected_het_positions.begin(), expected_het_positions.end());
    assert(tempSet == heterozygous_positions);
    assert_msg(true, "PhaseSetComputer", "Heterozygous positions.");

    delete superreads;
    delete accessible_positions;
}

void test_non_overlapping_reads() {
    ReadSet* superreads = mock_superreads();
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_3();
    std::vector<uint32_t> accessible_positions = {100, 200, 300, 400, 600, 700, 800, 900, 1100, 1200, 1300, 1400, 1500};

    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(true);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(true);
    
    {
        read1->addVariant(100, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(200, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        
        read2->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read2->addVariant(400, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read2->addVariant(500, std::vector<uint32_t>{}, Entry::BLANK);
        read2->addVariant(600, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
    }

    compute_phasesets(&accessible_positions, read_set, superreads);
    
    assert(read1->hasPhaseSet());
    assert(read1->getPhaseSet() == 100);
    assert(read2->hasPhaseSet());
    assert(read2->getPhaseSet() == 400);

    assert_msg(true, "PhaseSetComputer", "Reads without overlapping HET position.");
}

void test_overlapping_reads() {

    ReadSet* superreads = mock_superreads();
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_3();
    std::vector<uint32_t> accessible_positions = {100, 200, 300, 400, 600, 700, 800, 900, 1100, 1200, 1300, 1400, 1500};

    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(true);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(true);
    
    {
        read1->addVariant(100, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(200, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(400, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);

        read2->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read2->addVariant(400, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read2->addVariant(500, std::vector<uint32_t>{}, Entry::BLANK);
        read2->addVariant(600, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
    }

    compute_phasesets(&accessible_positions, read_set, superreads);
    
    assert(read1->hasPhaseSet());
    assert(read1->getPhaseSet() == 100);
    assert(read2->hasPhaseSet());
    assert(read2->getPhaseSet() == 100);

    assert_msg(true, "PhaseSetComputer", "Reads with overlapping HET position.");

}

void test_read_not_covering_het() {
    ReadSet* superreads = mock_superreads();
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_3();
    std::vector<uint32_t> accessible_positions = {100, 200, 300, 400, 600, 700, 800, 900, 1100, 1200, 1300, 1400, 1500};

    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(true);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(true);
    
    {
        read1->addVariant(200, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        
        read2->addVariant(500, std::vector<uint32_t>{}, Entry::BLANK);
    }

    compute_phasesets(&accessible_positions, read_set, superreads);
    
    assert(!read1->hasPhaseSet());
    assert(!read2->hasPhaseSet());
    
    assert_msg(true, "PhaseSetComputer", "Reads not covering any HET positions.");
}

void test_unselected_read() {
    ReadSet* superreads = mock_superreads();
    std::vector<variant_information_t> variant_info_table = mock_variant_info_table_3();
    std::vector<uint32_t> accessible_positions = {100, 200, 300, 400, 600, 700, 800, 900, 1100, 1200, 1300, 1400, 1500};

    ReadSet* read_set = new ReadSet();
    Read* read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(true);
    Read* read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(true);
    Read* read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(false);
    
    {
        read1->addVariant(100, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(200, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        
        read2->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read2->addVariant(400, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read2->addVariant(500, std::vector<uint32_t>{}, Entry::BLANK);
        read2->addVariant(600, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);

        read3->addVariant(100, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read3->addVariant(200, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read3->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read3->addVariant(400, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
    }

    compute_phasesets(&accessible_positions, read_set, superreads);
    // Spans two different phaseblocks
    assert(!read3->hasPhaseSet());
    delete read_set;

    read_set = new ReadSet();
    read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(true);
    read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(true);
    read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(false);

    {
        read1->addVariant(100, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(200, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(400, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        
        read2->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read2->addVariant(400, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read2->addVariant(500, std::vector<uint32_t>{}, Entry::BLANK);
        read2->addVariant(600, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);

        read3->addVariant(100, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read3->addVariant(200, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read3->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read3->addVariant(400, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
    }
    compute_phasesets(&accessible_positions, read_set, superreads);
    assert(read3->hasPhaseSet());
    assert(read3->getPhaseSet() ==  100);
    delete read_set;

    read_set = new ReadSet();
    read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(true);
    read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(true);
    read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(false);

    {
        read1->addVariant(100, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(200, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        
        read2->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read2->addVariant(400, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read2->addVariant(500, std::vector<uint32_t>{}, Entry::BLANK);
        read2->addVariant(600, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);

        read3->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read3->addVariant(400, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
    }
    compute_phasesets(&accessible_positions, read_set, superreads);
    assert(read3->hasPhaseSet());
    assert(read3->getPhaseSet() ==  400);
    delete read_set;

    read_set = new ReadSet();
    read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(true);
    read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(true);
    read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(false);

    {
        read1->addVariant(100, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(200, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        
        read2->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read2->addVariant(400, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read2->addVariant(500, std::vector<uint32_t>{}, Entry::BLANK);
        read2->addVariant(600, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);

        read3->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
    }
    compute_phasesets(&accessible_positions, read_set, superreads);
    assert(!read3->hasPhaseSet());

    delete read_set;

    read_set = new ReadSet();
    read1 = new Read("read1", 60, 0); read_set->add(read1); read1->setSelected(true);
    read2 = new Read("read2", 60, 0); read_set->add(read2); read2->setSelected(true);
    read3 = new Read("read3", 60, 0); read_set->add(read3); read3->setSelected(false);

    {
        read1->addVariant(100, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(200, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read1->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        
        read2->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read2->addVariant(400, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read2->addVariant(500, std::vector<uint32_t>{}, Entry::BLANK);
        read2->addVariant(600, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);

        read3->addVariant(200, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read3->addVariant(300, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read3->addVariant(500, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
        read3->addVariant(700, std::vector<uint32_t>{}, Entry::EQUAL_SCORES);
    }
    compute_phasesets(&accessible_positions, read_set, superreads);
    assert(!read3->hasPhaseSet());
    
    assert_msg(true, "PhaseSetComputer", "Reads not selected for phasing.");
    
}

void test_phasesetcomputer() {
    
    test_heterozygous_positions();
    test_non_overlapping_reads();
    test_overlapping_reads();
    test_read_not_covering_het();
    test_unselected_read();
}