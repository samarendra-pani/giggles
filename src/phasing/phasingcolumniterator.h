/*
Taken from whatshap (version 2.8)
Original filename: src/columniterator.h
*/

#ifndef PHASING_COLUMNITERATOR_H
#define PHASING_COLUMNITERATOR_H

#include <vector>
#include <memory>
#include <list>

#include "../entry.h"
#include "../readset.h"
#include "../genotypingalgorithm.h"

class PhasingColumnIterator {
public:
	PhasingColumnIterator(const ReadSet& set, const std::vector<GenotypingAlgorithm::variant_information_t>* variant_info_table, bool is_first_phasing_round);
	~PhasingColumnIterator();
	/** Returns the total number of columns, i.e. the number of columns
	 *  that will be returned by get_next. */
	uint32_t get_column_count(); 
	/** Returns the total number of reads. */
	uint32_t get_read_count(); 
	bool has_next();
	/** Ownership of Entry objects remains with the PhasingColumnIterator. Pointers
	 *  remain valid only until iterator is destructed. */
	std::unique_ptr<std::vector<const Entry*> > get_next();
	const std::vector<uint32_t>* get_positions();
	/** Moves iterator such that next call to get_next() will return 
	 *  column k. */
	void jump_to_column(size_t k);

private:
	typedef struct active_read_t {
		size_t read_index;
		size_t active_entry;
		active_read_t(size_t read_index) : read_index(read_index), active_entry(0) {}
		active_read_t(size_t read_index, size_t active_entry) : read_index(read_index), active_entry(active_entry) {}
	} active_read_t;
	
	const ReadSet& set;
	/** The number of columns already written. */
	size_t n;
	/** Index of the read that is to be examined next. */
	size_t next_read_index;
	std::list<active_read_t> active_reads;
	std::vector<Entry*> blank_entries;
	// positions of the variants
	std::vector<uint32_t>* positions;
	// number of active alleles at each position
	std::vector<uint32_t>* n_active_alleles;
	// if position is a structural variant
	std::vector<bool>* sv_flag;
	// first_reads[k] is the index of the first read (i.e. lowest index) active at column k,
	// in case no read is active in column k, then first_reads[k] is the index of the first read
	// that will become active after column k.
	std::vector<size_t> first_reads;
	// flag indicating if this is the first phasing round. Not considering SVs in the first round.
	bool is_first_phasing_round;
};

#endif
