/*
Taken from whatshap (version 2.8)
Original filename: src/pedigreedptable.h
*/

#ifndef PHASING_DP_TABLE_H
#define PHASING_DP_TABLE_H

#include <array>
#include <vector>
#include <memory>

#include "phasingcolumnindexingscheme.h"
#include "phasingcolumniterator.h"
#include "../entry.h"
#include "../read.h"
#include "../readset.h"
#include "../vector2d.h"
#include "../genotypelikelihoods.h"
#include "../genotypingalgorithm.h"

class GenotypingAlgorithm;

class PhasingDPTable {
	private:
		// pointer to the read set
		ReadSet* read_set;
		/* pointer to the variant information table from genotyping algorithm.
		* contains information about each variant position, including its genotype likelihoods.
		*/
		const std::vector<GenotypingAlgorithm::variant_information_t>* variant_info_table;
		// vector of indexingschemes
		std::vector<PhasingColumnIndexingScheme*> indexers;
		// optimal score and its index in the rightmost DP table column
		uint32_t optimal_score;
		uint32_t optimal_score_index;
		// projection_column_table[c] contains the projection column "between" columns c and c+1
		std::vector<std::vector<uint32_t>* > projection_column_table;
		// index_backtrace_table[c][i][t] indicates the index (=bipartition) in column c from which the
		// i-th entry in the FORWARD projection of column c comes from, assuming a transmission value of t
		std::vector<std::vector<uint32_t>* > index_backtrace_table;
		PhasingColumnIterator input_column_iterator;
		// optimal path obtained from backtrace
		std::vector<uint32_t> index_path;

		// helper function to pull read ids out of read column
		std::unique_ptr<std::vector<uint32_t> > extract_read_ids(const std::vector<const Entry *>& entries);

		/** Initializes/clears all member variables associated with the DP table, i.e. indexers, index_backtrace_table,
		 *  transmission_backtrace_table, optimal_score, optimal_score_index, optimal_transmission_value, and previous_transmission_value. */
		void clear_table();
		void compute_table();
		/** Computes the DP column at the given index, assuming that the previous column
		 *  has already been computed. */
		void compute_column(size_t column_index, std::unique_ptr<std::vector<const Entry*>> current_input_column = nullptr);

		/** Returns the number of set bits. */
		static size_t popcount(size_t x);

		template <class T>
		void init(std::vector<T*>& v, size_t size) {
			for(size_t i=0; i<v.size(); ++i) {
				if (v[i] != nullptr) {
					delete v[i];
				}
			}
			v.assign(size, nullptr);
		}

	public:
		/** Constructor.
		 *  @param read_set DP table is constructed for the contained reads. Ownership is retained
		 *                  by caller. Pointer must remain valid during the lifetime of this PhasingDPTable.
		 *  @param variant_info_table Contains information about each variant position, including its genotype likelihoods.
		 *  @param first_phasing_round Indicates whether this is the first phasing round (not considering SVs) or not (considering SVs).
		 *  
		 */
		PhasingDPTable(ReadSet* read_set, const std::vector<GenotypingAlgorithm::variant_information_t>* variant_info_table, bool first_phasing_round);
	
		~PhasingDPTable();

		uint32_t get_optimal_score();

		/** Computes optimal haplotypes and adds them (in the form of "super reads") to 
		 *  the given read_set.
		 *
		 *   @param output_read_set Must have as many entries as there are individuals. The haplotypes for individual
		 *                          with index i in the pedigree (given at construction time) are added to output_read_set->at(i).
		 */
		void get_super_reads(ReadSet* output_read_set);

		/** Performs a backtrace through the DP table and returns optimal partitioning of the reads.
		 *  Pointer ownership is transferred to caller. */
		std::vector<bool>* get_optimal_partitioning();
};

#endif
