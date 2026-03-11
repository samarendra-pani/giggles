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
void find_phasesets_tag_reads(ComponentFinder<uint32_t>* component_finder, const std::unordered_set<uint32_t>* accessible_positions_set, ReadSet* read_set, const std::unordered_set<uint32_t>& heterozygous_positions);

/**
 * Original python function name in whatshap/cli/phase.py: compute_overall_components
 */

/**
 * Finds the heterozygous positions from the superreads and calls find_phasesets_tag_reads.
 * Variants are considered to be in the same component if a read exists that covers both.
 * A component is identified by the representative position (usually the smallest/leftmost variant).
 * @param accessible_positions List of all variant positions that were phased by the DP table (have 2 active alleles).
 * @param read_set Pointer to the set of reads containing variant information. All the reads are in this object.
 * @param superreads Pointer to the set of superreads (should contain exactly two reads).
 * 
 * @return Nothing. The reads in read_set are tagged with their phaseset ID.
 */
void compute_phasesets(const std::vector<uint32_t>* accessible_positions, ReadSet* read_set, ReadSet* superreads);

#endif