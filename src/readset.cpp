// Code taken from WhatsHap (https://github.com/whatshap/whatshap)

#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <unordered_set>
#include <iomanip>

#include "readset.h"

using namespace std;

ReadSet::ReadSet() {
}


ReadSet::~ReadSet() {
	for (size_t i=0; i<reads.size(); ++i) {
		delete reads[i];
	}
}


void ReadSet::add(Read* read) {
	name_and_source_id_t name_and_source_id = name_and_source_id_t(read->getName(), read->getSourceID());
	if (read_name_map.find(name_and_source_id) != read_name_map.end()) {
		throw std::runtime_error("ReadSet::add: duplicate read name.");
	}
	reads.push_back(read);
	read_name_map[name_and_source_id] = reads.size() - 1;
}


string ReadSet::toString() {
	ostringstream oss;
	oss << "ReadSet:" << endl;
	for (size_t i=0; i<reads.size(); ++i) {
		oss << "  " << setw(5) << i << ' ' << reads[i]->toString() << endl;
	}
	return oss.str();
}


void ReadSet::initialize() {
	// Sort the reads by position
	std::sort(reads.begin(), reads.end(), read_comparator_t());
	// Update read_name_map
	read_name_map.clear();
	pos_to_entry_map.clear();
	for (size_t i=0; i<reads.size(); ++i) {
		Read* read = reads[i];
		read_name_map[name_and_source_id_t(read->getName(), read->getSourceID())] = i;
		read->setID(i);
		// iterating through the Entries and creating the map from variant position to Entry pointers.
		for (size_t j=0; j < read->getVariantCount(); ++j) {
			uint32_t pos = read->getPosition(j);
            Entry* e = read->getEntry(j);
            auto [it, inserted] = pos_to_entry_map.try_emplace(pos, std::vector<Entry*>());
            it->second.push_back(e);
		}
	}
}


vector<uint32_t>* ReadSet::get_positions() const {
	unordered_set<uint32_t> position_set;
	for (size_t i=0; i<reads.size(); ++i) {
		reads[i]->addPositionsToSet(&position_set);
	}
	vector<uint32_t>* positions = new vector<uint32_t>(position_set.begin(), position_set.end());
	std::sort(positions->begin(), positions->end());
	return positions;
}


uint32_t ReadSet::size() const {
	return reads.size();
}


Read* ReadSet::get(uint32_t i) const {
	return reads[i];
}


Read* ReadSet::getByName(std::string name, int source_id) const {
	read_name_map_t::const_iterator it = read_name_map.find(name_and_source_id_t(name,source_id));
	if (it == read_name_map.end()) {
		return 0;
	} else {
		return reads[it->second];
	}
}


ReadSet* ReadSet::subset(const IndexSet* indices) const {
	ReadSet* result = new ReadSet();
	IndexSet::const_iterator it = indices->begin();
	for (; it != indices->end(); ++it) {
		result->add(new Read(*(reads[*it])));
	}
	return result;
}


void ReadSet::assign_selection_status(const IndexSet* indices) {
	IndexSet::const_iterator it = indices->begin();
	std::unordered_set<int> index_set;
	for (; it != indices->end(); ++it) {
		reads[*it]->setSelected(true);
	}
}


void ReadSet::resetTags() {
	for (size_t i = 0; i < reads.size(); ++i) {
		reads[i]->resetTags();
	}
}

void ReadSet::setEntryAlleles(uint32_t pos, std::vector<bool> active_alleles) {
	for (Entry* e: pos_to_entry_map[pos]) {
		e->set_allele_type(active_alleles);
	}
}

std::vector<Entry*> ReadSet::TEST_get_pos_to_entry_map(uint32_t pos) {
	return pos_to_entry_map[pos];
}