#ifndef GENOTYPE_HMM
#define GENOTYPE_HMM

#include <array>
#include <vector>
#include <memory>

#include "haplotypemapper.h"
#include "column.h"
#include "columniterator.h"
#include "entry.h"
#include "read.h"
#include "readset.h"
#include "vector2d.h"
#include "backwardcolumniterator.h"
#include "transitionprobabilitycomputer.h"
#include "emissionprobabilitycomputer.h"
#include "variantinfo.h"

class GenotypeHMM {
	private:
		// number of haplotypes present in the graph.
		uint32_t num_haplotypes;
		
		/**
		 * contains information about each variant position which is to be genotyped.
		 * The information includes the position, reference allele, and alternate alleles.
		 */ 
		std::vector<variant_information_t>* variant_info_table;
		
		// the input sequencing reads
		ReadSet* read_set;
		
		/**
		 * the recombination cost vector based on Li Stephens model
		 * @see UniformRecombinationCostComputer in giggles/giggles/utils.py
		 */
		const std::vector<float>& recombcost;
		
		// storing the Columns that are made from read clusters at variant positions
		std::vector<Column*> hmm_columns;

		/**
		 * maps between haplotype pairs to the reduced states selected based on previous genotyping
		 */
		std::vector<HaplotypeMapper*> haplotype_mapper_table;
		
		/**
		 * tables storing the values calculated from the backward pass
		 * table[idx] stores the value at column idx.
		 * 
		 * Note: we do not need one for forward pass since we do not store
		 * those values. We create forward pass values and directly calculate
		 * genotype likelihoods based on stored backward pass.
		 */
		std::vector<std::vector<long double>* > backward_pass_table;
		
		/**
		 * vector to store the forward probabilities for the current and previous positions.
		 */
		std::vector<long double> previous_forward_probabilities;
		std::vector<long double> current_forward_probabilities;

		/**
		 * vectors for storing the helper variables of previous column
		 */
		std::vector<long double> alpha_helper_1;	// This helper value is the alpha(*,*) value.
		std::vector<std::vector<long double>> alpha_helper_2;	// This helper value is the alpha(R1,*) value.
		std::vector<std::vector<long double>> alpha_helper_3;	// This helper value is the alpha(*,R2) value.
		/**
		 * vectors for storing the helper variables of current column
		 */
		std::vector<long double> curr_alpha_helper_1;
		std::vector<std::vector<long double>> curr_alpha_helper_2;
		std::vector<std::vector<long double>> curr_alpha_helper_3;

		//iterator used to iterate the columns of the input matrix (forward)
		ColumnIterator input_column_iterator;
		
		// iterator used to iterate the columns of the input matrix (backward)
		BackwardColumnIterator backward_input_column_iterator;
		
		// scaling parameters
		std::vector<long double> scaling_parameters;

		// helper to pull read ids out of read column
		std::unique_ptr<std::vector<uint32_t> > extract_read_ids(const std::vector<const Entry *>& entries);
		
		/**
		 * clears the backward pass tables used for the HMM
		 */
		void clear_backward_table();
		
		/**
		 * Forward Pass
		 * computes the forward probabilities using compute_forward_column()
		 */
		void compute_forward_prob();
		
		/**
		 * Backward Pass
		 * computes the backward probabilities using compute_backward_column()
		 */
		void compute_backward_prob();
		
		/**
		 * Computes various data structures that are needed for the Genotyping HMM.
		 * - The Column which store the active reads/read clusters.
		 * - The HaplotypeMapper which determines what haplotype pair states are active.
		 */
		void compute_index();

		/**
		 * Computes the forward probabilities for column_index
		 * 
		 * This function iterates over the states in column_index and calculates the forward probabilities for
		 * these states using the forward probabilities from column_index - 1.
		 * 
		 * We calculate the finally genotype likelihoods in this function (since backward probabilities have been
		 * calculated prior to executing this).
		 */
		void compute_forward_column(size_t column_index, std::unique_ptr<std::vector<const Entry*>> current_input_column = nullptr);

		/**
		 * Computes the backward probabilities for column_index - 1 (NOTE: for column_index - 1 and not column_index)
		 * 
		 * This function iterates over the states in column_index and calculates the contribution of these
		 * states to the states in column_index - 1
		 * 
		 * While during the iteration of column_index, the backward values for column_index - 1 are calculated,
		 * they are NOT normalized.
		 * When the function is called for column_index - 1, the values are normalized (and are now probabilities).
		 */
		void compute_backward_column(size_t column_index, std::unique_ptr<std::vector<const Entry*>> current_input_column = nullptr);

		/**
		 * Given bipartition index (b_index), the index inside the bipartition (r_index),
		 * and number of states within each bipartition (num_states),
		 * returns the index of the state being referred to.
		 */
		uint32_t get_node_index(uint32_t b_index, uint32_t r_index, uint32_t num_states);

		// used to initialize/clear tables
		template<class T>
		void init(std::vector<T*>& v, size_t size)
		{
			for(size_t i=0; i<v.size(); ++i) {
					if (v[i] != nullptr) {
							delete v[i];
					}
			}
			v.assign(size,nullptr);
		}

	public:
		/** Constructor
		 * @param read_set   DP table is constructed for the given reads. Ownership is retained by caller.
		 *			Pointer must remain valid during the lifetime of this GenotypeDPTable.
		 * @param recombcost phred scaled recombination probabilities
		 * @param n_references number of reference haplotypes in the graph
		 * @param variant_info_table contains information about each variant position which is to be genotyped.
		 */
		GenotypeHMM(ReadSet* read_set, const std::vector<float>& recombcost, const uint32_t& n_references, std::vector<variant_information_t>* variant_info_table);
		~GenotypeHMM();

		// returns the computed genotype likelihoods for a given position
		std::vector<long double> get_genotype_likelihoods(uint32_t position);

};
#endif
