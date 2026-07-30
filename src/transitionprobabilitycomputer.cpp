#include "transitionprobabilitycomputer.h"
#include <cassert>

TransitionProbabilities calculate_transition_probabilities(
    uint32_t varpos1,
    uint32_t varpos2,
    long double transition_constant,
    uint32_t num_haplotypes) {

    const long double s = static_cast<long double>(num_haplotypes);
    assert(varpos2 > varpos1);
    const long double r = ((long double)(varpos2 - varpos1))*transition_constant;
    
    const long double exponent = expl(-(r / s));
    const long double p = (1.0L - exponent) / s;
    const long double q = exponent + p;
    
    return { 
        p * p, // p2
        q * q, // q2
        p * q  // pq
    };
}