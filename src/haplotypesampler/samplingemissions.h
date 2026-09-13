/**
 * Original implementation at https://github.com/eblerjana/pangenie/blob/master/src/samplingemissions.hpp (Commit f682fb6)
 */

#ifndef SAMPLING_EMISSIONS_H
#define SAMPLING_EMISSIONS_H

#include "../entry.h"
#include <memory>
#include <vector>

class SamplingEmissions {
public:
	SamplingEmissions(const std::vector<const Entry*>& entries, uint32_t n_alleles);
	unsigned int get_emission_cost(unsigned int allele_id) const;
	void penalize(unsigned int allele_id, unsigned int penalty);
private:
	std::vector<unsigned int> allele_penalties;
	unsigned int default_penalty;

};

#endif // SAMPLING_EMISSIONS