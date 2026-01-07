/*
* Functions taken from WhatsHap (version 2.8)
* Original filename: whatshap/cli/phase.py
* Code adapted to C++ and modified.
*/

#ifndef PHASESET_COMPUTER_H
#define PHASESET_COMPUTER_H

#include <algorithm>
#include <cassert>

#include "componentfinder.h"
#include "../../readset.h"

/**
 * Original python function name in whatshap/cli/phase.py: find_components
 */

/**
 * Finds connected components of variants based on read coverage.
 * * It tags the reads in the readset with their phaseset ID.
 * Variants are considered to be in the same component if a read exists that covers both.
 * A component is identified by the representative position (usually the smallest/leftmost variant).
 * Phaseblocks are determined by reads selected for phasing (i.e., those that contributed to the DP table).
 * The unselected reads are then tagged with the phaseset if all their variants belong to the same component.
 * @param component_finder Pointer to the ComponentFinder object used to manage components.
 * @param accessible_positions_set Set of all variant positions that were phased by the DP table.
 * @param read_set Pointer to the set of reads containing variant information. All the reads are in this object.
 * @param heterozygous_positions List of positions to restrict component building. 
 * If empty, all variants in reads are used. 
 * If not empty, only variants at these positions are used to link components.
 * * @return Nothing. The reads in read_set are tagged with their phaseset ID.
 */
void find_phasesets_tag_reads(ComponentFinder<uint32_t>* component_finder, const std::unordered_set<uint32_t>* accessible_positions_set, ReadSet* read_set, const std::unordered_set<uint32_t>& heterozygous_positions) {
    
    
    bool filter_by_het = !heterozygous_positions.empty();
    
    for (uint32_t i = 0; i < read_set->size(); ++i) {
        Read* read = read_set->get(i);
        if (!read->isSelected()) {
            continue; // Skip unselected reads.
        }
        std::vector<uint32_t> read_positions;
        read_positions.reserve(read->getVariantCount());
        for (uint32_t j = 0; j < read->getVariantCount(); ++j) {
            uint32_t pos = read->getPosition(j);

            // Check if it is a phased position
            if (accessible_positions_set->find(pos) == accessible_positions_set->end()) { continue; }
            // Check heterozygous constraint
            // If filter_by_het is false, we skip this check
            if (filter_by_het && heterozygous_positions.find(pos) == heterozygous_positions.end()) { continue; }
            read_positions.push_back(pos);
        }

        // Merge components for all pairs of positions in the read
        if (read_positions.size() > 1) {
            uint32_t first = read_positions[0];
            for (size_t m = 1; m < read_positions.size(); ++m) {
                component_finder->merge(first, read_positions[m]);
            }
        }
    }

    // tagging reads with their phaseset ID (representative position)
    for (uint32_t i = 0; i < read_set->size(); ++i) {
        Read* read = read_set->get(i);
        if (read->isSelected()) {
            uint32_t ps = component_finder->find(read->firstPosition());
            read->setPhaseSet(ps);
            continue;
        }
        else {
            // make sure all positions in the read belong to the same component
            uint32_t rep = component_finder->find(read->firstPosition());
            bool all_same_component = true;
            for (uint32_t j = 1; j < read->getVariantCount(); ++j) {
                if (rep != component_finder->find(read->getPosition(j))) {
                    all_same_component = false;
                    break;
                }
            }
            if (all_same_component) {
                read->setPhaseSet(rep);
            }
        }  
    }
}

/**
 * Original python function name in whatshap/cli/phase.py: compute_overall_components
 */

/**
 * Finds the heterozygous positions from the superreads and calls find_phasesets_tag_reads.
 * Variants are considered to be in the same component if a read exists that covers both.
 * A component is identified by the representative position (usually the smallest/leftmost variant).
 * @param accessible_positions List of all variant positions that were phased by the DP table.
 * @param read_set Pointer to the set of reads containing variant information. All the reads are in this object.
 * @param superreads Pointer to the set of superreads (should contain exactly two reads).
 * * @return Nothing. The reads in read_set are tagged with their phaseset ID.
 */
void compute_phasesets(const std::vector<uint32_t>* accessible_positions, ReadSet* read_set, ReadSet* superreads) {
    
    std::unordered_set<uint32_t> heterozygous_positions;
    std::unordered_set<uint32_t> accessible_positions_set(accessible_positions->begin(), accessible_positions->end());

    assert (superreads->size() == 2); // two superreads represeting the two haplotypes
    Read* superread0 = superreads->get(0);
    Read* superread1 = superreads->get(1);

    assert (superread0->getVariantCount() == superread1->getVariantCount()); // both superreads should have same number of variants

    for (uint32_t i = 0; i < superread0->getVariantCount(); ++i) {
        assert (superread0->getPosition(i) == superread1->getPosition(i)); // both superreads should have variants at same positions
        // skip positions that are not accessible
        if (accessible_positions_set.find(superread0->getPosition(i)) == accessible_positions_set.end()) { continue; }
        
        Entry::allele_t allele0 = superread0->getEntry(i)->get_allele_type();
        Entry::allele_t allele1 = superread1->getEntry(i)->get_allele_type();

        if ((allele0 == Entry::ALLELE1 && allele1 == Entry::ALLELE2) ||
            (allele0 == Entry::ALLELE2 && allele1 == Entry::ALLELE1)) {
            heterozygous_positions.insert(superread0->getPosition(i));
        }
    }
    ComponentFinder<uint32_t> component_finder(*accessible_positions);
    find_phasesets_tag_reads(&component_finder, &accessible_positions_set, read_set, heterozygous_positions);
}

#endif