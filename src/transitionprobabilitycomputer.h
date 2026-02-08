#ifndef TRANSITIONPROBABILITYCOMPUTER_H
#define TRANSITIONPROBABILITYCOMPUTER_H

/**
 * Calculates the pr and qr values of the Li-Stephens transition probility model.
 * pr -> result.first
 * qr -> result.second
 * 
 * Note:
 * We do not normalize these values since we do an overall normalization at the end.
 */
std::pair<long double, long double>  calculate_transition_probabilities(const float& recombcost, const uint32_t& num_haplotypes) {
    long double s = (long double)num_haplotypes;
    long double r = (long double)recombcost;
    
    std::pair<long double, long double> result;
    result.first = (1.0 - exp(-(r/s))) / s;         // calculating pr
    result.second = exp(-(r/s)) + result.first;     // calculating qr
}

#endif // TRANSITIONPROBABILITYCOMPUTER_H
