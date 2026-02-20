#ifndef TRANSITIONPROBABILITYCOMPUTER_H
#define TRANSITIONPROBABILITYCOMPUTER_H

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

#endif // TRANSITIONPROBABILITYCOMPUTER_H
