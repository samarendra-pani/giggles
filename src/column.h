#ifndef COLUMN_H
#define COLUMN_H

#include <vector>
#include <memory>
#include "columnindexingiterator.h"
#include "readset.h"

class ColumnIndexingIterator;

class Column {
private:
	// read IDs of active reads at the variant position
	std::vector<unsigned int> read_ids;
	// read IDS of untagged active reads at the variant position
	std::vector<unsigned int> untagged_read_ids;
	// read IDs of active reads at the next variant position
	std::vector<unsigned int> next_read_ids;
	// read IDS of untagged active reads at the next variant position
	std::vector<unsigned int> untagged_next_read_ids;
	// read IDs of active, non terminating reads at the variant position
	std::vector<unsigned int> act_nonterminating_read_ids;
	// read IDs of active, terminating reads at the variant position
	std::vector<unsigned int> act_terminating_read_ids;
	unsigned int n_references;
	
public:

	Column(const unsigned int index, const unsigned int* n_ref, const std::vector<unsigned int>& read_ids, const std::vector<unsigned int>& next_read_ids, ReadSet* set);
	
	// returns the column size
	unsigned int get_column_size();

	// return a pointer to the read ids
	std::vector<unsigned int> * get_read_ids();

	// return a pointer to the active non terminating read ids of the column
	std::vector<unsigned int> * get_active_nonterminating_read_ids();

	// return a pointer to the active terminating read ids of the column
	std::vector<unsigned int> * get_active_terminating_read_ids();

	// return a pointer to the active non terminating read ids of the column
	std::vector<unsigned int> * get_next_read_ids();

	// returns a pointer to the bipartition iterator which uses graycode.
	std::unique_ptr<ColumnIndexingIterator> get_iterator(ReadSet* set);

	// gets the index value using the bipartition index and reference index;
	unsigned int get_index(unsigned int b_index, unsigned int r_index);

	// returns the compatible bipartitions of bipartition b_index (of pos v+1) in column v
	std::vector<unsigned int> get_backward_compatible_bipartitions(int b_index, ReadSet* set);

};

#endif
