// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#ifndef READ_H
#define READ_H

#include <string>
#include <vector>
#include <unordered_set>

#include "entry.h"

class Read {
public:
	Read(const std::string& name, int mapq, int source_id, int sample_id, int reference_start = -1, const std::string& BX_tag = "", int reg_const = 10, double base_const = 2.718);
	virtual ~Read() {}
	std::string toString();
	void addHaplotag(std::string hp, int ps);
    int getHaplotag() const;
    int getPhaseSet() const;
	bool hasHaplotag() const;
	bool hasPhaseSet() const;
	void addVariant(int position, int allele, std::vector<double> em, int quality);
	void sortVariants();
	/** Returns the position of the first variant. **/
	int firstPosition() const;
	/** Returns the position of the last variant. **/
	int lastPosition() const;
	void setID(int id);
	int getID() const;
	/** Add all positions contained in this read to the given set. */
	void addPositionsToSet(std::unordered_set<unsigned int>* set);
	int getPosition(size_t variant_idx) const;
	void setPosition(size_t variant_idx, int position);
	int getAllele(size_t variant_idx) const;
	void setAllele(size_t variant_idx, int allele);
	std::vector<long double> getEmissionProbability(size_t variant_idx) const;
	void setEmissionProbability(size_t variant_idx, std::vector<double> emission);
	int getQuality(size_t variant_idx) const;
	void setQuality(size_t variant_idx, int quality);
	const Entry* getEntry(size_t variant_idx) const;
	int getVariantCount() const;
	const std::string& getName() const;
	const std::vector<int>& getMapqs() const;
	void addMapq(int mapq);
	int getSourceID() const;
	int getSampleID() const;
	int getReferenceStart() const;
	const std::string& getBXTag() const;
	int getRegConst() const;
	double getBaseConst() const;
	bool isSorted() const;
	bool hasBXTag() const;
	
	
private:
	typedef struct enriched_entry_t {
		Entry entry;
		int position;
		enriched_entry_t(int position, int allele, std::vector<double> em, int quality, int reg_const, double base_const) :
			entry(0,allele,em,quality, reg_const, base_const), position(position) {}
	} enriched_entry_t;
	
	typedef struct entry_comparator_t {
		entry_comparator_t() {}
		bool operator()(const enriched_entry_t& e1, const enriched_entry_t& e2) {
			return e1.position < e2.position;
		}
	} entry_comparator_t;

	std::string name;
	std::vector<int> mapqs;
	int source_id;
	int sample_id;
	int id;
	int reference_start;
	std::string BX_tag;
	std::vector<enriched_entry_t> variants;
	int hp;
	int ps;
	int reg_const;
	double base_const;
};

#endif
