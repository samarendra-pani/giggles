// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#ifndef READ_H
#define READ_H

#include <string>
#include <vector>
#include <unordered_set>

#include "entry.h"

class Read {
public:
	Read(const std::string& name, uint32_t mapq, uint32_t source_id, int reference_start = -1);
	virtual ~Read() {}
	std::string toString();
	
	void addVariant(uint32_t position, uint32_t allele, std::vector<uint32_t> scores);
	void addHaplotag(std::string hp, int ps);
    /** Add all positions contained in this read to the given set. */
	void addPositionsToSet(std::unordered_set<uint32_t>* set);
	void addMapq(uint32_t mapq);
	
	const std::string& getName() const;
	uint32_t getSourceID() const;
	uint32_t getID() const;
	const Entry* getEntry(size_t variant_idx) const;
	uint32_t getPosition(size_t variant_idx) const;
	std::vector<uint32_t> getScores(size_t variant_idx) const;
	uint32_t getAllele(size_t variant_idx) const;
	int getHaplotag() const;
    int getPhaseSet() const;
	int getReferenceStart() const;
	const std::vector<uint32_t>& getMapqs() const;
	uint32_t getVariantCount() const;
	/** Returns the position of the first variant. **/
	uint32_t firstPosition() const;
	/** Returns the position of the last variant. **/
	uint32_t lastPosition() const;
	
	void setID(uint32_t id);
	void setPosition(size_t variant_idx, uint32_t position);
	void setScores(size_t variant_idx, std::vector<uint32_t> scores);
	void setAllele(size_t variant_idx, uint32_t allele);
	

	bool hasHaplotag() const;
	bool hasPhaseSet() const;
	
	void sortVariants();
	bool isSorted() const;
	
	bool isSelected() const; // is this read selected for phasing?
	void setSelected(bool selected); // set whether this read is selected for phasing
	
	
private:
	typedef struct enriched_entry_t {
		uint32_t position; // position on the reference
		uint32_t index; // zero-based index for multi-allelic variants
		Entry entry; // the record entry
		enriched_entry_t(uint32_t position, uint32_t allele, std::vector<uint32_t> scores) :
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
