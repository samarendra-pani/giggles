#ifndef HAPLOTAG_COMPUTER_H
#define HAPLOTAG_COMPUTER_H

#include <cassert>

#include "../readset.h"

/**
 * Function to calculate the edit distance between a read and a superread (haplotype) inferered by the DP table.
 * * Since some positions in the read may not have been phased (due to being an SV), only positions present in the superread are considered.
 * @param read Pointer to the read whose distance from the superread is to be calculated.
 * @param superread Pointer to the superread representing a haplotype.
 * @param position_to_index Mapping from variant position to index in the superread.
 * @return The edit distance between the read and the superread.
 */
uint32_t calculate_distance_from_superread(Read* read, Read* superread, const std::unordered_map<uint32_t, uint32_t>& position_to_index) {

    uint32_t distance = 0;
    for (uint32_t i = 0; i < read->getVariantCount(); ++i) {
        uint32_t pos = read->getPosition(i);
        auto it = position_to_index.find(pos);
        if (it != position_to_index.end()) {
            uint32_t superread_index = it->second;
            Entry* read_entry = read->getEntry(i);
            Entry* superread_entry = superread->getEntry(superread_index);
            if (read_entry->get_allele_type() != superread_entry->get_allele_type()) {
                distance += 1; // Increment distance for mismatch
            }
        }
    }
    return distance;
}

//TODO: check if superreads actually correspong to H1 and H2
/**
 * Haplotag unselected reads based on their distance to the two superreads representing the haplotypes.
 * @param read_set Pointer to the set of reads to be haplotagged.
 * @param superreads Pointer to the set of superreads (should contain exactly two reads).
 * @return Nothing. The reads in read_set are haplotagged.
 */
void haplotag_unselected_reads(ReadSet* read_set, ReadSet* superreads) {
    assert (superreads->size() == 2); // two superreads represeting the two haplotypes
    Read* superread0 = superreads->get(0);
    Read* superread1 = superreads->get(1);

    assert (superread0->getVariantCount() == superread1->getVariantCount()); // both superreads should have same number of variants

    // Create a mapping from position to index for superreads
    std::unordered_map<uint32_t, uint32_t> position_to_index;
    for (uint32_t i = 0; i < superread0->getVariantCount(); ++i) {
        uint32_t pos = superread0->getPosition(i);
        position_to_index[pos] = i;
    }

    // Haplotag each read based on distance to superreads
    for (uint32_t i = 0; i < read_set->size(); ++i) {
        Read* read = read_set->get(i);
        if (read->isSelected()) {
            continue; // Skip selected reads. Their haplotag comes from the DP table.
        }
        if (read->getPhaseSet() == -1) {
            continue; // Skip reads without a phaseset.
        }
        uint32_t distance_to_hap0 = calculate_distance_from_superread(read, superread0, position_to_index);
        uint32_t distance_to_hap1 = calculate_distance_from_superread(read, superread1, position_to_index);

        if (distance_to_hap0 < distance_to_hap1) {
            read->addHaplotag("H1");
        } else if (distance_to_hap1 < distance_to_hap0) {
            read->addHaplotag("H2");
        }
    }
}

/**
 * Haplotag reads that were selected for the phasing algorithm, based on the given partitioning.
 * @param read_set Pointer to the set of reads to be haplotagged.
 * @param partitioning Pointer to a vector of booleans indicating the partitioning of the reads.
 * @return Nothing. The reads in read_set are haplotagged.
 */
void haplotag_selected_reads(ReadSet* read_set, const std::vector<bool>* partitioning) {
    assert (read_set->size() == partitioning->size());
    for (uint32_t i = 0; i < read_set->size(); ++i) {
        Read* read = read_set->get(i);
        if (!read->isSelected()) {
            assert (partitioning->at(i) == false); // they should be partitioned in the DP table as false/
            continue;
        }
        if (partitioning->at(i)) {
            read->addHaplotag("H1");
        } else {
            read->addHaplotag("H2");
        }
    }
}

#endif