// Code taken from WhatsHap (https://github.com/whatshap/whatshap)

#ifndef COLUMNITERATOR_H
#define COLUMNITERATOR_H

#include <vector>
#include <memory>
#include <list>

#include "entry.h"
#include "readset.h"
#include "variantinfo.h"

/**
 * 
 */
class ColumnIterator {
	public:
		ColumnIterator(const ReadSet& set, const std::vector<variant_information_t>* variant_info_table);
		~ColumnIterator();
		/** 
		 * Returns the total number of columns, i.e. the number of columns
		 * that will be returned by get_next.
		 */
		uint32_t get_column_count(); 
		/** Returns the total number of reads. */
		uint32_t get_read_count(); 
		// Checks if more columns are available in the forward context
		bool has_next();
		// Checks if more columns are available in the backward context
		bool has_prev();
		/** 
		 * Returns a pointer to the active column in forward context.
		 */
		std::unique_ptr<std::vector<const Entry*> > get_next();
		/** 
		 * Returns a pointer to the active column in backward context.
		 */
		std::unique_ptr<std::vector<const Entry*> > get_prev();
		/** 
		 * Moves iterator such that next call to get_next()/get_prev() will return 
		 *  column k.
		 * 
		 * Sets both n and m to k.
		 */
		void jump_to_column(uint32_t k);

	private:
		typedef struct active_read_t {
			size_t read_index;
			size_t active_entry;
			active_read_t(size_t read_index) : read_index(read_index), active_entry(0) {}
			active_read_t(size_t read_index, size_t active_entry) : read_index(read_index), active_entry(active_entry) {}
		} active_read_t;
		
		const ReadSet& set;
		// The number of columns already written in forward context.
		uint32_t n;
		// The number of columns already written in backward context.
		uint32_t m;
		// Index of the read that is to be examined next.
		size_t next_read_index;
		std::list<active_read_t> active_reads;
		// Blank entries for the current column
		std::vector<Entry*> current_blank_entries;
		const std::vector<variant_information_t>* variant_info_table;
		/**
		 * first_reads[k] is the index of the first read (i.e. lowest index) active at column k,
		 * in case no read is active in column k, then first_reads[k] is the index of the first read
		 * that will become active after column k.
		 */
		std::vector<size_t> first_reads;
		/**
		 * Gets the active reads at column k
		 */
		void get_active_reads(uint32_t k);
};

#endif
