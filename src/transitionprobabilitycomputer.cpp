#include "transitionprobabilitycomputer.h"
#include <cmath>
#include <iostream>
#include <cassert>

#include "phredgenotypelikelihoods.h"

using namespace std;

TransitionProbabilityComputer::TransitionProbabilityComputer(const float& recombcost, const vector<int>& next_allele_reference) {
    int s = next_allele_reference.size();
    float r = recombcost;
    this->pr = (1.0 - exp(-(r/s))) / s;
    this->qr = exp(-(r/s)) + pr;
}
