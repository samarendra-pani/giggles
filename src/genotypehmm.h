#ifndef GENOTYPE_HMM
#define GENOTYPE_HMM

#include <array>
#include <vector>
#include <memory>

#include "column.h"
#include "columniterator.h"
#include "entry.h"
#include "read.h"
#include "readset.h"
#include "vector2d.h"
#include "backwardcolumniterator.h"
#include "transitionprobabilitycomputer.h"
#include "genotypingalgorithm.h"

class GenotypingAlgorithm;

class GenotypeHMM {
	private:
		// number of reference samples
		unsigned int n_references;
		
		/* variant information table
		* contains information about each variant position which is to be genotyped.
		* The information includes the position, reference allele, and alternate alleles.
		*/ 
		std::vector<GenotypingAlgorithm::variant_information_t>* variant_info_table;
		
		// the input sequencing reads
		ReadSet* read_set;
		
		// the recombination cost vector
		const std::vector<float>& recombcost;
		
		// indexing schemes
		std::vector<Column*> hmm_columns;
		
		// projection_column_table[c] contains the projection column between columns c and c+1
		std::vector<std::vector<long double>* > forward_pass_column_table;
		std::vector<std::vector<long double>* > backward_pass_column_table;

		std::vector<std::vector<unsigned int>* > active_reads;
		
		//iterator used to iterate the columns of the input matrix (forward)
		ColumnIterator input_column_iterator;
		
		// iterator used to iterate the columns of the input matrix (backward)
		BackwardColumnIterator backward_input_column_iterator;
		
		// stores the transmission probability computers for each column. object at index i contains the probability computer between index i and i+1.
		std::vector<TransitionProbabilityComputer*> transition_probability_table;
		
		// scaling parameters
		std::vector<long double> scaling_parameters;

		// helper to pull read ids out of read column
		std::unique_ptr<std::vector<unsigned int> > extract_read_ids(const std::vector<const Entry *>& entries);
		
		// initializes all members associated with the DP table
		void clear_forward_table();
		void clear_backward_table();
		
		// forward pass: computes the forward probabilities
		void compute_forward_prob();
		
		// backward pass: computes the backward probabilities
		void compute_backward_prob();
		
		// computes the index for each column
		void compute_index();

		// computes column of forward probabilities of given index, assuming previous column was already computed (from left to right)
		void compute_forward_column(size_t column_index, std::unique_ptr<std::vector<const Entry*>> current_input_column = nullptr);

		// computes column of backward probabilities of given index, assuming previous column was already computed (from right to left)
		void compute_backward_column(size_t column_index, std::unique_ptr<std::vector<const Entry*>> current_input_column = nullptr);

		// returns the number of bits set
		static size_t popcount(size_t x);

		// updates the emission probabilities
		void update_emission_probability(Vector2D<long double>* em_prob, const int bit_changed, const ColumnIndexingIterator& iterator, std::vector<const Entry *>& entries);

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
		GenotypeHMM(ReadSet* read_set, const std::vector<float>& recombcost, const unsigned int& n_references, std::vector<GenotypingAlgorithm::variant_information_t>* variant_info_table);
		~GenotypeHMM();

		// returns the computed genotype likelihoods for a given position
		std::vector<long double> get_genotype_likelihoods(unsigned int position);

};
#endif
