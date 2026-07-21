#include "haplotagcomputer.h"

uint32_t calculate_distance_from_superread(Read* read, Read* superread, const std::unordered_map<uint32_t, uint32_t>& position_to_index) {

    uint32_t distance = 0;
    for (uint32_t i = 0; i < read->getVariantCount(); ++i) {
        uint32_t pos = read->getPosition(i);
        auto it = position_to_index.find(pos);
        if (it != position_to_index.end()) {
            uint32_t superread_index = it->second;
            Entry* read_entry = read->getEntry(i);
            Entry* superread_entry = superread->getEntry(superread_index);
            /** Not considering positions where the allele type is BLANK or EQUAL SCORE */
            if (superread_entry->get_allele_type() == Entry::BLANK || superread_entry->get_allele_type() == Entry::EQUAL_SCORES) {
                continue;
            }
            if (read_entry->get_allele_type() == Entry::BLANK || read_entry->get_allele_type() == Entry::EQUAL_SCORES) {
                continue;
            }
            assert(read_entry->get_allele_type() == Entry::ALLELE1 || read_entry->get_allele_type() == Entry::ALLELE2);
            assert(superread_entry->get_allele_type() == Entry::ALLELE1 || superread_entry->get_allele_type() == Entry::ALLELE2);
            if (read_entry->get_allele_type() != superread_entry->get_allele_type()) {
                distance += 1; // Increment distance for mismatch
            }
        }
    }
    return distance;
}

void haplotag_unselected_reads(ReadSet* read_set, ReadSet* superreads) {
    assert (superreads->size() == 2); // two superreads represeting the two haplotypes
    Read* superread0 = superreads->get(0);
    Read* superread1 = superreads->get(1);

    assert (superread0->getVariantCount() == superread1->getVariantCount()); // both superreads should have same number of variants

    // Create a mapping from position to index for superreads
    std::unordered_map<uint32_t, uint32_t> position_to_index;
    for (uint32_t i = 0; i < superread0->getVariantCount(); ++i) {
        uint32_t pos = superread0->getPosition(i);
        assert(pos == superread1->getPosition(i));
        position_to_index[pos] = i;
    }

    // Haplotag each read based on distance to superreads
    size_t count = 0;
    for (uint32_t i = 0; i < read_set->size(); ++i) {
        Read* read = read_set->get(i);
        /** Skip selected reads. Their haplotag comes from the DP table. */
        if (read->isSelected()) { continue; }
        /** Skip reads without a phaseset. */
        if (!read->hasPhaseSet()) { continue; }
        uint32_t distance_to_hap0 = calculate_distance_from_superread(read, superread0, position_to_index);
        uint32_t distance_to_hap1 = calculate_distance_from_superread(read, superread1, position_to_index);
        if (distance_to_hap0 < distance_to_hap1) {
            count++;
            read->setHaplotag(false);
        } else if (distance_to_hap1 < distance_to_hap0) {
            count++;
            read->setHaplotag(true);
        }
    }
    std::cerr << "[Core::Phasing] Haplotagged " << count << " unselected reads (out of " << read_set->size() << " total reads)" << std::endl;
}

void haplotag_selected_reads(ReadSet* read_set, const std::vector<bool>* partitioning) {
    assert (read_set->size() == partitioning->size());
    size_t count = 0;
    for (uint32_t i = 0; i < read_set->size(); ++i) {
        Read* read = read_set->get(i);
        if (!read->isSelected()) {
            assert (partitioning->at(i) == false); // they should be partitioned in the DP table as false
            continue;
        }
        count++;
        read->setHaplotag(partitioning->at(i));
    }
    std::cerr << "[Core::Phasing] Haplotagged " << count << " selected reads (out of " << read_set->size() << " total reads)" << std::endl;
}