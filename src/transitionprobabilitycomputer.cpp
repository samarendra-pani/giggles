#include "transitionprobabilitycomputer.h"

TransitionProbabilities calculate_transition_probabilities(float recombcost, uint32_t num_haplotypes) {
    // Using 'L' suffix for long double literals to maintain precision
    const long double s = static_cast<long double>(num_haplotypes);
    const long double r = static_cast<long double>(recombcost);
    
    const long double exponent = expl(-(r / s));
    const long double p = (1.0L - exponent) / s;
    const long double q = exponent + p;

    // Direct return for Mandatory Copy Elision
    return { 
        p * p, // p2
        q * q, // q2
        p * q  // pq
    };
}