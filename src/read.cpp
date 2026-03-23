// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <cassert>


#include "read.h"

using namespace std;

Read::Read(const std::string& name, uint32_t mapq, uint32_t source_id): 
	name(name),
	mapqs(1, mapq),
	source_id(source_id) {
	this->id = -1;
	selected = false;
	hp = false;
	has_hp = false;
	ps = 0;
	has_ps = false;
	is_clustered = false;
	cluster_id = 0;
	constrained_cluster_id = 0;
	has_constrained_cluster = false;
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


void Read::setHaplotag(std::string hp) {
	if (hp == "H1") {this->hp = false; this->has_hp = true;}
	if (hp == "H2") {this->hp = true; this->has_hp = true;}
	//if (hp == "none") {throw std::runtime_error("Read with 'none' haplotag found. These should be filtered.");}
	if (hp == "none") {this->hp = false; this->has_hp = false;}
}

void Read::setClusterID(uint32_t cluster_id) {
	this->cluster_id = cluster_id;
}

void Read::setClusterStatus(bool is_clustered) {
	this->is_clustered = is_clustered;
}

bool Read::isClustered() const {
	return is_clustered;
}

uint32_t Read::getClusterID() const {
	return cluster_id;
}

bool Read::getClusterStatus() const {
	return is_clustered;
}

void Read::setConstrainedClusterID(uint32_t constrained_cluster_id) {
	this->has_constrained_cluster = true;
	this->constrained_cluster_id = constrained_cluster_id;
}

uint32_t Read::getConstrainedClusterID() const {
	if (!has_constrained_cluster) {
		return (uint32_t)-1;
	}
	return constrained_cluster_id;
}

bool Read::hasConstrainedCluster() const {
	return has_constrained_cluster;
}

void Read::unsetPhaseSet() {
	this->ps = 0;
	this->has_ps = false;
}

void Read::setPhaseSet(uint32_t ps) {
	this->has_ps = true;
	this->ps = ps;
}

bool Read::getHaplotag() const {
	if (has_hp == false) {
		throw std::runtime_error("Haplotag not set for read " + name);
	}
	return hp;
}

uint32_t Read::getPhaseSet() const {
	if (has_ps == false) {
		throw std::runtime_error("Phase set not set for read " + name);
	}
	return ps;
}

bool Read::hasHaplotag() const {
	return has_hp;
}

bool Read::hasPhaseSet() const {
	return has_ps;
}

void Read::addVariant(uint32_t position, vector<uint32_t> scores) {
	variants.push_back(enriched_entry_t(position, scores));
}

void Read::addVariant(uint32_t position, vector<long double> scores) {
	variants.push_back(enriched_entry_t(position, scores));
}

void Read::addVariant(uint32_t position, vector<uint32_t> scores, Entry::allele_t allele) {
	variants.push_back(enriched_entry_t(position, scores, allele));
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


std::vector<long double> Read::getEmissionScores(size_t variant_idx) const {
	assert(variant_idx < variants.size());
	return variants[variant_idx].entry.get_emission_scores();
}


void Read::setScores(size_t variant_idx, std::vector<uint32_t> scores) {
	assert(variant_idx < variants.size());
	variants[variant_idx].entry.set_scores(scores);
}

void Read::setEmissionScores(size_t variant_idx, std::vector<long double> scores) {
	assert(variant_idx < variants.size());
	variants[variant_idx].entry.set_emission_scores(scores);
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