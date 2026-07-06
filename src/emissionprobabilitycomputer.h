#ifndef EMISSIONPROBABILITYCOMPUTER_H
#define EMISSIONPROBABILITYCOMPUTER_H

#include "vector2d.h"
#include "bipartitioniterator.h"
#include "entry.h"

/**
 * @brief Computes the emission probability matrix
 *
 * Given some bipartition of reads at a variant position, this matrix stores at position i, j
 * the probability of the reads in biparition 1 having allele i and reads in biparition 2 having allele j.
 */
class EmissionProbabilityComputer {

    public:

        EmissionProbabilityComputer(uint32_t n_alleles);
        
        /**
         * Get the emission probability for allele pair (i, j)
         */
        long double at(uint32_t i, uint32_t j) const;

        /**
         * Given a bit that has changed in the bipartition defined by the BiparitionIterator,
         * we update the emission probabilities.
         * What value does emission_probability_table(i, j) have?
         * 
         * INITIALIZATION: if (bit_changed < 0) --> this indicates that no bit was changed and that the bipartion is the first bipartition we are looking at.
         *      - We look at the bit assignment of all the clusters which are in .
         *      - The bipartition is defined by the b_index.
         *      - Based on the bipartion of the clusters, we find the reads inside the clusters and assign then emission score i if they are bipartition 0 otherwise emission score j.
         * UPDATE: if (bit_changed >= 0) --> this indicates that some bit was changed.
         *      - Since Initialization already created a starting emission probability table, this code block just updates the table.
         *      - The bit that was flippped corresponds to a cluster. We extract the reads that were present in the cluster.
         *      - If read was flipped from bipartition 0 to 1, we divide the table value with emission score i and multiply with emission score j
         */
        void update_emission_probability(const int bit_changed, const BipartitionIterator& iterator, std::vector<const Entry *>& entries);


    private:

        /**
         * A 2D table of size n x n where n is the number of alleles at that position.
         * The value at (i, j) gives the following value:
         * Given a bipartition of reads B = {B0, B1} where B0 corresponds to reads assigned to haplotype 0 and B1 corresponds to reads assigned to haplotype 1.
         * We have the assignment of reads in B0 to allele i and reads in B1 to allele j.
         * The value at (i, j) gives the probability of observing the reads in their assigned alleles given this bipartition and allele assignment.
         */
        Vector2D<long double> emission_probability_table;

        /**
         * An unordered map which stores map between index of reads (in the set of active reads at the position) and their new bipartition.
         * 
         * This is passed to BipartitionIterator as a reference where this is calculated.
         */
        std::unordered_map<uint32_t, bool> changed_reads;

        /**
         * Helper variables to speed up emission update.
         */
        std::vector<long double> i_multiplier;
        std::vector<long double> j_multiplier;

};


#endif // EMISSIONPROBABILITYCOMPUTER_H