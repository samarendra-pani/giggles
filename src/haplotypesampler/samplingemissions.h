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
	unsigned int get_emission_cost(unsigned short allele_id) const;
	void penalize(unsigned short allele_id, unsigned short penalty);
private:
	std::vector<unsigned short> allele_penalties;
	unsigned int default_penalty;

};

#endif // SAMPLING_EMISSIONS