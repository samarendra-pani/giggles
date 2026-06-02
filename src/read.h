// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#ifndef READ_H
#define READ_H

#include <string>
#include <vector>
#include <unordered_set>
#include <cstdint>

#include "entry.h"

class Read {
	
	public:
		Read(const std::string& name, uint32_t mapq, uint32_t source_id);
		virtual ~Read() {}
		std::string toString();

		// adding a variant
		void addVariant(uint32_t position, std::vector<uint32_t> scores);
		void addVariant(uint32_t position, std::vector<long double> scores);
		// adding a variant whose allele has been pre-computed. used for the superread haplotypes.
		void addVariant(uint32_t position, std::vector<uint32_t> scores, Entry::allele_t allele);
		void setHaplotag(bool hp);
		void unsetHaplotag();
		void unsetPhaseSet();
		void setPhaseSet(uint32_t ps);
		void setClusterID(uint32_t cluster_id);
		void setClusterStatus(bool is_clustered);
		void setConstrainedClusterID(uint32_t constrained_cluster_id);
		/** Add all positions contained in this read to the given set. */
		void addPositionsToSet(std::unordered_set<uint32_t>* set);
		void addMapq(uint32_t mapq);

		const std::string& getName() const;
		uint32_t getSourceID() const;
		uint32_t getID() const;
		Entry* getEntry(size_t variant_idx);
		uint32_t getPosition(size_t variant_idx) const;
		std::vector<long double> getEmissionScores(size_t variant_idx) const;
		bool getHaplotag() const;
		uint32_t getPhaseSet() const;
		uint32_t getClusterID() const;
		uint32_t getConstrainedClusterID() const;
		const std::vector<uint32_t>& getMapqs() const;
		uint32_t getVariantCount() const;
		/** Returns the position of the first variant. **/
		uint32_t firstPosition() const;
		/** Returns the position of the last variant. **/
		uint32_t lastPosition() const;

		void setID(uint32_t id);
		void setPosition(size_t variant_idx, uint32_t position);
		void setScores(size_t variant_idx, std::vector<uint32_t> scores);
		void setEmissionScores(size_t variant_idx, std::vector<long double> scores);


		bool hasHaplotag() const;
		bool hasPhaseSet() const;
		bool isClustered() const;
		bool hasConstrainedCluster() const;

		void sortVariants();
		bool isSorted() const;

		bool isSelected() const; // is this read selected for phasing?
		void setSelected(bool selected); // set whether this read is selected for phasing


	private:
		typedef struct enriched_entry_t {
			uint32_t position; // position on the reference
			uint32_t index; // zero-based index for multi-allelic variants
			Entry entry; // the record entry
			enriched_entry_t(uint32_t position, std::vector<uint32_t> scores) :
				entry(0, scores), position(position) {}
			enriched_entry_t(uint32_t position, std::vector<long double> scores) :
				entry(0, {}), position(position) { entry.set_emission_scores(scores); }
			enriched_entry_t(uint32_t position, std::vector<uint32_t> scores, Entry::allele_t allele) :
				entry(0, scores), position(position) { entry.set_allele_type(allele); }
		} enriched_entry_t;

		typedef struct entry_comparator_t {
			entry_comparator_t() {}
			bool operator()(const enriched_entry_t& e1, const enriched_entry_t& e2) {
				if (e1.position != e2.position) {
					return e1.position < e2.position;
				}
				return e1.index < e2.index;
			}
		} entry_comparator_t;

		std::string name;
		std::vector<uint32_t> mapqs;
		uint32_t source_id;
		uint32_t id;
		std::vector<enriched_entry_t> variants;
		
		bool selected; // selected for phasing
		
		/**
		 * clustering information based on phasing
		 */
		bool has_hp; // stores whether the read has a haplotag
		bool hp; // haplotag value - false = 0, true = 1
		bool has_ps; // stores whether the read has a phaseset
		uint32_t ps;  // phaseset value
		bool is_clustered; // whether the read has been assigned to a cluster
		uint32_t cluster_id; // cluster ID of the read
		uint32_t constrained_cluster_id; // cluster ID which cannot be in the same bipartition as the cluster_id
		bool has_constrained_cluster; // whether the read has a constrained cluster
};

#endif
