// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#include <cassert>
#include "columnindexingiterator.h"
#include "column.h"
#include <math.h>

using namespace std;

Column::Column(const unsigned int index, const unsigned int* n_ref, const std::vector<unsigned int>& read_ids, const std::vector<unsigned int>& next_read_ids, ReadSet* set) : read_ids(read_ids), next_read_ids(next_read_ids) {

	// Finding the untagged read ids from the variant position
	for (auto read = begin(read_ids); read != end(read_ids); read++) {
		if (!set->get(*read)->hasHaplotag()) {
			untagged_read_ids.push_back(*read);
		}
	}
	// Finding the untagged read ids from the next variant position
	for (auto read = begin(next_read_ids); read != end(next_read_ids); read++) {
		if (!set->get(*read)->hasHaplotag()) {
			untagged_next_read_ids.push_back(*read);
		}
	}
	bool read_found;
	n_references = *n_ref; // each reference sample has 2 haplotype paths.
	for (auto read = begin(read_ids); read != end(read_ids); read++) {
		read_found = false;
		for (auto next_read = begin(next_read_ids); next_read != end(next_read_ids); next_read++){
			if (*read == *next_read) {
				act_nonterminating_read_ids.push_back(*read);
				read_found = true;
				break;
			}
		}
		if (!read_found) {
			act_terminating_read_ids.push_back(*read);
		}
	}
}

unsigned int Column::get_index(unsigned int b_index, unsigned int r_index) {
	unsigned int index = (n_references*n_references*b_index)+r_index;
	assert (index < this->get_column_size());
	return index;
}

unique_ptr<ColumnIndexingIterator> Column::get_iterator(ReadSet* set) {
	return unique_ptr<ColumnIndexingIterator>(new ColumnIndexingIterator(this, set));
}

unsigned int Column::get_column_size() {
	return pow(2 ,untagged_read_ids.size()) * pow(n_references,2);
}

vector<unsigned int> * Column::get_read_ids() {
	return &(this->read_ids);
}

vector <unsigned int> * Column::get_active_nonterminating_read_ids() {
	return &(this->act_nonterminating_read_ids);
}

vector <unsigned int> * Column::get_active_terminating_read_ids() {
	return &(this->act_terminating_read_ids);
}

vector <unsigned int> * Column::get_next_read_ids() {
	return &(this->next_read_ids);
}

// Gives the compatible bipartitions in the left column (at position index) given the bipartition index of a bipartition of the right column (at position index + 1)
vector<unsigned int> Column::get_backward_compatible_bipartitions(int b_index, ReadSet* set) {
	vector<unsigned int> compatible_bipartition = {0};
	int base = 0;
    int count = 0;
	// This now works with the untagged read IDs (and not with all the reads)
	// This for loop deals with all the reads in read_ids that are less than max(next_read_ids)
	// Example: Let read_ids = {1,2,3,4,5,6,7,8,9,10} and next_read_ids = {2,3,5,6,8}
	// So this for loop creates the necessary bipartitions by processing read_ids till Read 8.
	for (int i = 0; i < untagged_next_read_ids.size(); i++) {
        int ri = untagged_next_read_ids.at(i);
		// This break statement if the last count++ step happened outside the while loop
		if (count >= untagged_read_ids.size()) break;
        while (ri > untagged_read_ids.at(count)) {
            if (set->get(untagged_read_ids.at(count))->getHaplotag() == -1) {
				compatible_bipartition.resize(2*compatible_bipartition.size());
				for (int j = 0; j < compatible_bipartition.size()/2; j++) {
					compatible_bipartition.at((compatible_bipartition.size()/2)+j) = compatible_bipartition.at(j) + pow(2, count);
				}
			}
			else {
				base += pow(2, count)*set->get(untagged_read_ids.at(count))->getHaplotag();
			}		
            count++;
			// Break out of while loop (otherwise read_ids.at(count) is not defined)
			if (count >= untagged_read_ids.size()) break;
        }
		// This break statement if the last count++ step happened inside the while loop and we don't want the base value to increase.
		if (count >= untagged_read_ids.size()) break;
        base = base + (pow(2,count)*(b_index%2));
        b_index = b_index/2;
        count++;
    }
	// Continuing the example from above.
	// The rest of the reads in read_ids which are 9 and 10 are processed here and these reads automatically lead to an increase in number of compatible bipartitions.
	for (int i = count; i < untagged_read_ids.size(); i++) {
		if (set->get(untagged_read_ids.at(i))->getHaplotag() == -1) {
			compatible_bipartition.resize(2*compatible_bipartition.size());
			for (int j = 0; j < compatible_bipartition.size()/2; j++) {
				compatible_bipartition.at((compatible_bipartition.size()/2)+j) = compatible_bipartition.at(j) + pow(2, i);
			}
		}
		else {
			base += pow(2, i)*set->get(untagged_read_ids.at(i))->getHaplotag();
		}
	}
	
	for (int i = 0; i < compatible_bipartition.size(); i++) {
        compatible_bipartition.at(i) += base;
    }
	return compatible_bipartition;
}
