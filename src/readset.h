// Code taken from WhatsHap (https://github.com/whatshap/whatshap)

#ifndef READSET_H
#define READSET_H

#include <string>
#include <vector>
#include <unordered_map>

#include "read.h"
#include "indexset.h"

class ColumnIterator;

class ReadSet {
public:
	ReadSet();
	virtual ~ReadSet();
	/** Ownership of pointer is transferred from caller to the ReadSet. */
	void add(Read* read);
	/** 
	 * - Sort reads by first variant position.
	 * - Assigns read_ids to all instances of Entry stored in the reads such that
	 * 	  each read_id matches the index of the corresponding read in the ReadSet.
	 * - Create the position to Entry map
	 */
	void initialize();
	/** Returns the set of SNP positions. To create this set,
	 *  this method iterates over all contained reads.
	 *  Caller owns the returned pointer. */
	std::vector<uint32_t>* get_positions() const;
	uint32_t size() const;
	std::string toString();
	/** Access a read in the set. Ownership stays with the ReadSet. */
	Read* get(uint32_t i) const;
	/** Access a read in the set by its name. Ownership stays with the ReadSet. */
	Read* getByName(std::string name, int source_id) const;
	/** Creates a subset of reads as given by the set of indices. Note that this
	 *  creates a COPY of each read.
	 */
	ReadSet* subset(const IndexSet* indices) const;
	/** Marks reads as selected/unselected based on the given indices. */
	void assign_selection_status(const IndexSet* indices);
	/** Resets all haplotags and cluster information */
	void resetTags();
	/* Sets the allele type for Entry objects at phasable positions */
	void setEntryAlleles(uint32_t pos, std::vector<bool> active_alleles);
	/* TEST function to inspect pos_to_entry_map */
	std::vector<Entry*> TEST_get_pos_to_entry_map(uint32_t pos);


private:
	typedef struct read_comparator_t {
		read_comparator_t() {}
		bool operator()(const Read* r1, const Read* r2) {
			// 1. Sort by presence of variants (reads with 0 variants go first)
			if (r1->getVariantCount() == 0 && r2->getVariantCount() > 0) return true;
			if (r2->getVariantCount() == 0 && r1->getVariantCount() > 0) return false;

			// 2. Standard case: sort by positions (if both have variants)
			if (r1->getVariantCount() > 0 && r2->getVariantCount() > 0) {
				if (r1->firstPosition() != r2->firstPosition()) {
					return r1->firstPosition() < r2->firstPosition();
				}
			}

			// 3. Break ties using string names directly (Deterministic & easy to test!)
			int name_cmp = r1->getName().compare(r2->getName());
			if (name_cmp != 0) {
				return name_cmp < 0; // Negative means r1 < r2
			}

			// 4. Ultimate tie-breaker
			return r1->getSourceID() < r2->getSourceID();
		}
	} read_comparator_t;

	typedef struct name_and_source_id_t {
		name_and_source_id_t(std::string name, int source_id) : name(name), source_id(source_id) {}
		bool operator==(const name_and_source_id_t& other) const {
			return (name.compare(other.name) == 0) && (source_id == other.source_id);
		}
		std::string name;
		int source_id;
	} name_and_source_id_t;

	typedef struct name_and_source_id_hasher_t {
		std::size_t operator()(const name_and_source_id_t& x) const {
			return (std::hash<std::string>()(x.name)) ^ (std::hash<int>()(x.source_id));
		}
	} name_and_source_id_hasher_t;

	std::vector<Read*> reads;
	// Maps names of reads it their index in the "reads" vector
	typedef std::unordered_map<name_and_source_id_t,size_t,name_and_source_id_hasher_t> read_name_map_t;
	read_name_map_t read_name_map;
	typedef std::unordered_map<uint32_t, std::vector<Entry*>> entry_pos_map_t;
	entry_pos_map_t pos_to_entry_map;
};

#endif