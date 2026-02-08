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
TransitionProbabilities calculate_transition_probabilities(const float& recombcost, const uint32_t& num_haplotypes) {
    long double s = (long double)num_haplotypes;
    long double r = (long double)recombcost;
    
    long double p = (1.0 - exp(-(r/s))) / s;   // calculating pr
    long double q = exp(-(r/s)) + p;           // calculating qr

    TransitionProbabilities result;
    result.p2 = pow(p, 2);
    result.q2 = pow(q, 2);
    result.pq = p*q;
}

#endif // TRANSITIONPROBABILITYCOMPUTER_H
