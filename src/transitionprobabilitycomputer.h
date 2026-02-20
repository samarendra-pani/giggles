#ifndef TRANSITIONPROBABILITYCOMPUTER_H
#define TRANSITIONPROBABILITYCOMPUTER_H

#include <cmath>
#include <cstdint>

struct TransitionProbabilities {
    long double p2;
    long double q2;
    long double pq;
};

/**
 * Calculates the pr and qr values of the Li-Stephens transition probility model: pr and qr.
 * Then it stores it in the struct TransitionProbabilities.
 *  
 * Note:
 * We do not normalize these values since we do an overall normalization at the end.
 */
TransitionProbabilities calculate_transition_probabilities(float recombcost, uint32_t num_haplotypes);

#endif // TRANSITIONPROBABILITYCOMPUTER_H
