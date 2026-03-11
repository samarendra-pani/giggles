#ifndef HAPLOTAG_COMPUTER_H
#define HAPLOTAG_COMPUTER_H

#include <cassert>

#include "../../readset.h"

/**
 * Function to calculate the edit distance between a read and a superread (haplotype) inferered by the DP table.
 * * Since some positions in the read may not have been phased (due to being an SV), only positions present in the superread are considered.
 * @param read Pointer to the read whose distance from the superread is to be calculated.
 * @param superread Pointer to the superread representing a haplotype.
 * @param position_to_index Mapping from variant position to index in the superread.
 * @return The edit distance between the read and the superread.
 */
uint32_t calculate_distance_from_superread(Read* read, Read* superread, const std::unordered_map<uint32_t, uint32_t>& position_to_index);

//TODO: check if superreads actually correspong to H1 and H2
/**
 * Haplotag unselected reads based on their distance to the two superreads representing the haplotypes.
 * @param read_set Pointer to the set of reads to be haplotagged.
 * @param superreads Pointer to the set of superreads (should contain exactly two reads).
 * @return Nothing. The reads in read_set are haplotagged.
 */
void haplotag_unselected_reads(ReadSet* read_set, ReadSet* superreads);

/**
 * Haplotag reads that were selected for the phasing algorithm, based on the given partitioning.
 * @param read_set Pointer to the set of reads to be haplotagged.
 * @param partitioning Pointer to a vector of booleans indicating the partitioning of the reads.
 * @return Nothing. The reads in read_set are haplotagged.
 */
void haplotag_selected_reads(ReadSet* read_set, const std::vector<bool>* partitioning);

#endif