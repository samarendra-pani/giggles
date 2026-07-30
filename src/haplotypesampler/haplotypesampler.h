/**
 * Original implementation at https://github.com/eblerjana/pangenie/blob/master/src/haplotypesampler.hpp (Commit f682fb6)
 */

#ifndef HAPLOTYPE_SAMPLER_H
#define HAPLOTYPE_SAMPLER_H

#include <vector>
#include <memory>
#include <stdexcept>
#include <cassert>
#include "samplingemissions.h"
#include "../variantinfo.h"
#include "../readset.h"


struct DPColumn {
	std::vector<uint32_t> column;
};

struct SampledPaths {
	std::vector<std::vector<uint32_t>> sampled_paths;
	/**
	* Given a column index, and a vector, mark indexes
	* occuring in this column as False.  
	**/
	std::vector<bool> mask_indexes(uint32_t column_index, uint32_t max_index) {
		std::vector<bool> masked(max_index+1, true);
		for (uint32_t i = 0; i < sampled_paths.size(); ++i) {
			if (column_index >= sampled_paths[i].size()) {
				throw std::runtime_error("HaplotypeSampler::SampledPaths::mask_indexes: column_index exceeds number of columns.");
			}
			uint32_t index = sampled_paths[i][column_index];
			if (index > max_index) {
				throw std::runtime_error("HaplotypeSampler::SampledPaths::mask_indexes: observed index exceeds max_index.");
			}
			masked[index] = false;
		}
		return masked;
	}

	/**
	* Given a column index and a path id, determine if there
	* was a recombination event between positions column_index-1 
	* and column_index.
	**/	
	bool recombination(uint32_t column_index, uint32_t path_id) {
		if (path_id >= sampled_paths.size()) {
			throw std::runtime_error("HaplotypeSampler::SampledPaths::recombination: path_id does not exist.");
		}

		if (column_index >= sampled_paths[path_id].size()) {
			throw std::runtime_error("HaplotypeSampler::SampledPaths::recombination: column_id does not exist.");
		}

		if (column_index > 0) {
			// recombination event if selected path changes
			return sampled_paths[path_id][column_index - 1] != sampled_paths[path_id][column_index];
		} else {
			// for first column, always set to false
			return false;
		}
	}
};


class HaplotypeSampler {
	
	public:
		/**
		* @param size size of the subsampled graph. Viterbi will be run this number of times
		* @param recombrate recombination rate
		* @param effective_N effective population size
		* @param best_scores vector in which DP score of each iteration is stored (mainly used for testing purposes)
		* @param remove_reference not consider reference sequence as an additional sampled path
		* @param path_output output paths of sampled path_ids to file
		* @param chromosome name of the chromosome (only used when writing path_output)
		* @param allele_penalty penality to penalize already covered alleles
		**/
		HaplotypeSampler(
			ReadSet* read_set,
			std::vector<variant_information_t>* variant_info_table,
			const uint32_t size, 
			const float recombrate = 1.26, 
			const float effective_N = 25000.0, 
			std::vector<uint32_t>* best_scores = nullptr, 
			const bool remove_reference = false, 
			const uint32_t allele_penalty = 10
		);

		// keeping it public for testing purposes ..
		void get_column_minima(std::vector<uint32_t>& column, std::vector<bool>& mask, uint32_t& first_id, uint32_t& second_id, uint32_t& first_val, uint32_t& second_val) const;

		// updates the variant information table to only contain the sampled paths.
		std::vector<variant_information_t> get_updated_variant_table(uint32_t ploidy);

		
	private:	

		/**
		* Compute emission probabilities from ReadSet. 
		*/
		void compute_emission_probabilities();

		/** Do one Viterbi pass and store the paths that have been used. 
		* This function can be applied several times and keeps track of used Viterbi path,
		* so that these nodes are ignored in the next pass (= call to this function). Used Viterbi
		* paths are stored in member variable "sampled_paths".
		**/
		void compute_viterbi_path(std::vector<uint32_t>* best_scores = nullptr);
		
		/**
		* Compute one column. Make sure that when iterating through the paths (= states), those paths
		* that were already used in previous Viterbi runs, are ignored. For this purpose, determine a list
		* of respective path ids first, so that these can be ignored while iterating the states.
		**/
		void compute_viterbi_column(uint32_t column_index);

		ReadSet* read_set;
		std::vector<variant_information_t>* variant_info_table;

		std::vector<DPColumn*> viterbi_columns;
		SampledPaths sampled_paths;
		std::vector< std::vector<uint32_t>* > viterbi_backtrace_columns;
		std::vector<bool> prev_mask;
		std::vector<SamplingEmissions> emission_costs;
		double recombrate;
		long double effective_N;
		uint32_t allele_penalty;
		uint32_t num_haplotypes;

		template<class T>
		void init(std::vector< T* >& c, uint32_t size) {
			for (uint32_t i = 0; i < c.size(); ++i) {
				if (c[i] != nullptr) delete c[i];
			}
			c.assign(size, nullptr);
		}
};

#endif // HAPLOTYPE_SAMPLER_HPP




/**
* Implementation ideas:
* - in each iteration of Viterbi, maxima need to be computed
* - this can be done efficiently by precomputing maxima for all sets: {p1, ..., pi-1, pi+1, ..., pn} (excluding each pi): 
*   compute overall maximum: x = max(p1,...,pn). Compute y = max(p1,...,px-1,px+1,...,pn). precomputed max for all sets with px is x, for all others y    
* - for each cell i compute: max(t0 * pi, t1*max(p1,...,pn)), the max(p1,..pn) term is precomputed
*
* Need: one class/function computing the transition probability / score /penalty, and one for the Emission.
* Transitions: compute phred-scaled Li stephans / WH transition score
* Emissions / local cost: 
**/