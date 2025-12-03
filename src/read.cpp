// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <cassert>


#include "read.h"

using namespace std;

Read::Read(const std::string& name, uint32_t mapq, uint32_t source_id, int reference_start): 
	name(name),
	mapqs(1, mapq),
	source_id(source_id),
	reference_start(reference_start) {
	this->id = -1;
	hp = -1;
	ps = -1;
}


string Read::toString() {
	ostringstream oss;
	oss << name << " mapq:(";
	for (size_t i=0; i<mapqs.size(); ++i) {
		if (i>0) oss << ",";
		oss << mapqs[i];
	}
	oss << ") source:" << source_id << " (";
	for (size_t i=0; i<variants.size(); ++i) {
		if (i>0) oss << ";";
		oss << "[" << variants[i].position << "," << variants[i].entry << "]";
	}
	oss << ")";
	return oss.str();
}


void Read::addHaplotag(std::string hp) {
	if (hp == "H1") {this->hp = 0;}
	if (hp == "H2") {this->hp = 1;}
	//if (hp == "none") {throw std::runtime_error("Read with 'none' haplotag found. These should be filtered.");}
	if (hp == "none") {this->hp = -1;}
}

void Read::addPhaseSet(int ps) {
	this->ps = ps;
}

int Read::getHaplotag() const {
	return hp;
}

int Read::getPhaseSet() const {
	return ps;
}

bool Read::hasHaplotag() const {
	return hp != -1;
}

bool Read::hasPhaseSet() const {
	return ps != -1;
}

void Read::addVariant(uint32_t position, vector<uint32_t> scores) {
	variants.push_back(enriched_entry_t(position, scores));
}

void Read::addVariant(uint32_t position, vector<uint32_t> scores, Entry::allele_t allele, uint32_t idx1, uint32_t idx2) {
	variants.push_back(enriched_entry_t(position, scores, allele, idx1, idx2));
}


void Read::sortVariants() {
	sort(variants.begin(), variants.end(), entry_comparator_t());
	for (size_t i=1; i<variants.size(); ++i) {
		if (variants[i-1].position == variants[i].position) {
			ostringstream oss;
			oss << "Duplicate variant in read " << name << " at position " << variants[i].position;
			throw std::runtime_error(oss.str());
		}
	}
}


uint32_t Read::firstPosition() const {
	if (variants.size() == 0) throw std::runtime_error("No variants present");
	return variants[0].position;
}


uint32_t Read::lastPosition() const {
	if (variants.size() == 0) throw std::runtime_error("No variants present");
	return variants[variants.size()-1].position;
}


void Read::setID(uint32_t id) {
	this->id = id;
	for (size_t i=0; i<variants.size(); ++i) {
		variants[i].entry.set_read_id(id);
	}
}


uint32_t Read::getID() const {
	return id;
}


void Read::addPositionsToSet(std::unordered_set<uint32_t>* set) {
	assert(set != 0);
	for (size_t i=0; i<variants.size(); ++i) {
		set->insert(variants[i].position);
	}
}


uint32_t Read::getPosition(size_t variant_idx) const {
	assert(variant_idx < variants.size());
	return variants[variant_idx].position;
}


void Read::setPosition(size_t variant_idx, uint32_t position) {
	assert(variant_idx < variants.size());
	variants[variant_idx].position = position;
}


std::vector<uint32_t> Read::getScores(size_t variant_idx) const {
	assert(variant_idx < variants.size());
	return variants[variant_idx].entry.get_scores();
}


void Read::setScores(size_t variant_idx, std::vector<uint32_t> scores) {
	assert(variant_idx < variants.size());
	variants[variant_idx].entry.set_scores(scores);
}


Entry* Read::getEntry(size_t variant_idx) {
	return &(variants[variant_idx].entry);
}


uint32_t Read::getVariantCount() const {
	return variants.size();
}


const string& Read::getName() const {
	return name;
}


const vector<uint32_t>& Read::getMapqs() const {
	return mapqs;
}


void Read::addMapq(uint32_t mapq) {
	mapqs.push_back(mapq);
}


uint32_t Read::getSourceID() const {
	return source_id;
}


int Read::getReferenceStart() const {
	return reference_start;
}

bool Read::isSorted() const {
	entry_comparator_t comparator;
	for (size_t i=1; i<variants.size(); ++i) {
		if (!comparator(variants[i-1],variants[i])) {
			return false;
		}
	}
	return true;
}

bool Read::isSelected() const {
	return selected;
}

void Read::setSelected(bool selected) {
	this->selected = selected;
}