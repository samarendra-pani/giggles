// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#ifndef READ_H
#define READ_H

#include <string>
#include <vector>
#include <unordered_set>

#include "entry.h"

class Read {
public:
	Read(const std::string& name, int mapq, int source_id, int reference_start = -1);
	virtual ~Read() {}
	std::string toString();
	void addHaplotag(std::string hp, int ps);
    int getHaplotag() const;
    int getPhaseSet() const;
	bool hasHaplotag() const;
	bool hasPhaseSet() const;
	void addVariant(int position, int allele, std::vector<uint32_t> scores);
	void sortVariants();
	/** Returns the position of the first variant. **/
	uint32_t firstPosition() const;
	/** Returns the position of the last variant. **/
	uint32_t lastPosition() const;
	void setID(uint32_t id);
	int getID() const;
	/** Add all positions contained in this read to the given set. */
	void addPositionsToSet(std::unordered_set<uint32_t>* set);
	uint32_t getPosition(size_t variant_idx) const;
	void setPosition(size_t variant_idx, uint32_t position);
	uint32_t getAllele(size_t variant_idx) const;
	void setAllele(size_t variant_idx, uint32_t allele);
	std::vector<uint32_t> getScores(size_t variant_idx) const;
	void setScores(size_t variant_idx, std::vector<uint32_t> scores);
	const Entry* getEntry(size_t variant_idx) const;
	uint32_t getVariantCount() const;
	const std::string& getName() const;
	const std::vector<uint32_t>& getMapqs() const;
	void addMapq(uint32_t mapq);
	uint32_t getSourceID() const;
	int getReferenceStart() const;
	bool isSorted() const;
	bool isSelected() const; // is this read selected for phasing?
	void setSelected(bool selected); // set whether this read is selected for phasing
	
	
private:
	typedef struct enriched_entry_t {
		int position; // position on the reference
		int index; // zero-based index for multi-allelic variants
		Entry entry; // the record entry
		enriched_entry_t(int position, int allele, std::vector<uint32_t> scores) :
			entry(0,allele,scores), position(position) {}
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
	int reference_start;
	std::vector<enriched_entry_t> variants;
	bool selected;
	int hp;
	int ps;
};

#endif
