#ifndef HAPLOTYPEMAPPER_H
#define HAPLOTYPEMAPPER_H

#include "genotypingalgorithm.h"
#include "vector2d.h"
#include "genotype.h"

/**
 * @brief Tracks the subset of haplotype pairs to consider for genotypes.
 * 
 * After one round of genotyping, the subsequent rounds of genotyping can use the previously
 * calculated genotype likelihoods to avoid unlikely genotypes.
 * HaplotypeMapper uses the genotype likelihoods and selects genotypes which are likely.
 * Then it maps the all the haplotype pairs that have the selected genotypes. 
 * 
 * This will be used to reduce the NxN space that each biparition (defined by BiparitionIterator)
 * in the HMM needs, where N is the number of haplotypes.
 * 
 * Let the number of selected haplotype pairs be M where NxN > M. In this class, we will refer to this
 * as "linearized reduced-space states"
 * 
 * @see BipartitionIterator
 * @note This class is hardcoded for ploidy 2
 */
class HaplotypeMapper {

    public:

        HaplotypeMapper(const GenotypeLikelihoods& genotype_likelihoods, const std::vector<int>& allele_references);

        // returns the size of the linearized reduced-space states, i.e. M
        uint32_t get_num_states() const;

        // get the index in linearized reduced-space states from the haplotype pair
        int get_state_index(uint32_t i, uint32_t j) const;

        // get the haplotype pair from the index in the linearized reduced-space
        std::pair<u_int32_t, u_int32_t> get_haplotypes_indices(uint32_t i) const;



    private:

        /**
         * This maps the NxN haplotype pairs to index of the linearized reduced-space states.
         * The vector stores -1 for the haplotype pairs which were not selected.
         */
        Vector2D<int> haplotype_to_selected_map;
        
        /**
         * The "Reverse Map" stores the map of the linearized reduced-space states to 
         * the original NxN matrix of haplotype pairs
         */
         std::vector<std::pair<u_int32_t, u_int32_t>> reverse_map;

};



# endif // HAPLOTYPEMAPPER_H