// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#ifndef GENOTYPE_H
#define GENOTYPE_H

#include<vector>
#include<set>
#include<string>
#include<algorithm>
#include <cstdint>

/**
 * Representation of a genotype of arbitrary ploidy with multi-allelic variants. Genotypes are 
 * unordered multisets of alleles, which are stored in 64bit word, divied into 16 chunks of 4 bits
 * each. The highest 4 bits (left most ones) encode the ploidy, ranging from 0 (invalid genotype) to
 * a maximum of 15. The remaining 15 4-bit-blocks encode the alleles in ascending order, i.e. the
 * lowest 4 bits (right most ones) encode the highest allele. There is a maximum of 16 different
 * alleles, which can be stored by this method.
 *
 * Genotypes have a canonical index in the VCF-format:
 * 
 * https://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes
 * 
 * Given the ploidy, there is a bijective mapping between genotypes and non-negative integer
 * numbers. In the diploid, bi-allelic case, the index is equal to the number of alternative
 * allles (either 0, 1 or 2):
 * 
 * 0 -> 0/0, 1 -> 0/1, 2 -> 1/1
 * 
 * For the bi-allelic case, this can easily be generalized for higher ploidy, e.g. 4:
 *
 * 0 -> 0/0/0/0, 1 -> 0/0/0/1, 2 -> 0/0/1/1, 3 -> 0/1/1/1, 4 -> 1/1/1/1
 *
 * For multi-allelic variants, there are more possible indices, which are enumerated in such way, that
 * first we have all genotypes, which only use the first allele (i.e. 0/0 for diploid genotypes. Then,
 * all genotypes with the first two alleles follow, then all genotypes using the first three alleles,
 * and so on. This way, the index can be interpreted without knowing the highest allele, as this can
 * be derived from the number itself (the ploidy is mandatory, though!):
 *
 * 0 -> 0/0
 * 1 -> 0/1, 2 -> 1/1
 * 3 -> 0/2, 4 -> 1/2, 5 -> 2/2
 * 6 -> 0/3, 7 -> 1/3, 8 -> 2/3, 9 -> 3/3
 * 
 * 0 -> 0/0/0/0
 * 1 -> 0/0/0/1, 2 -> 0/0/1/1, 3 -> 0/1/1/1, 4 -> 1/1/1/1
 * 5 -> 0/0/0/2, 6 -> 0/0/1/2, 7 -> 0/1/1/2, 8 -> 1/1/1/2, 9 -> 0/0/2/2, 10 -> 0/1/2/2, 11 -> 1/1/2/2, 12 -> 0/2/2/2, 13 -> 1/2/2/2, 14 -> 2/2/2/2
 * etc.
 * 
 * Note: The ploidy is hard set to 2 in this implementation. The code is computationally unfeasible for higher ploidy.
 */
class Genotype {
	public:
		
		// maximum supported number of alleles = 2^15
		const static uint32_t MAX_ALLELES = 32768;

		// maximum supported ploidy = 2
		const static uint32_t MAX_PLOIDY = 2;

		// creates an empty genotype with no alleles.
		Genotype();

		// creates a genotype of given ploidy using the canonical index (see class description).
		Genotype(uint32_t index, uint32_t ploidy);

		// creates a genotype from a list of given alleles.
		Genotype(std::vector<uint32_t> alleles);
	
		// returns the genotype's alleles as a vector.
		std::vector<uint32_t> as_vector() const;

		// returns whether the genotype is empty (i.e. invalid).
		bool is_none() const;
	
		// returns the canonical index of the genotype (see class description).
		uint32_t get_index() const;

		// returns the ploidy.
		uint32_t get_ploidy() const;

		// set the ploidy
		void set_ploidy(const uint32_t ploidy);

		// returns the genotype as readable string.
		std::string toString() const;

		// returns whether the genotype is homozygous.
		bool is_homozygous() const;

		// returns whether the genotype has ploidy 2 and only alleles 0 and 1.
		bool is_diploid_and_biallelic() const;

		// operators
		friend bool operator== (const Genotype &g1, const Genotype &g2);
		friend bool operator!= (const Genotype &g1, const Genotype &g2);
		friend bool operator< (const Genotype &g1, const Genotype &g2);

	private:
		/**
		 * Bitstring for storage. Each allele is encoded in 15 bits.
		 * We support ploidy of 1 (for chrY) or 2. 
		 * 		- we use 30 bits for alleles.
		 * 		- we use 2 bits for ploidy (00 is considered as uninitialized genotype and 11 is not allowed).
		 */
		uint32_t gt;
	
		// general manipulation methods
		uint32_t get_position(const uint32_t pos) const;
		void set_position(const uint32_t pos, const uint32_t allele);
};

/**
 * Creates a sorted vector of alleles from a given canonical index and ploidy.
 */
std::vector<uint32_t> convert_index_to_alleles(uint32_t index, uint32_t ploidy);

/**
 * Finds the cannonical index of the genotype using the sorted vector of alleles
 */
uint32_t convert_alleles_to_index(std::vector<uint32_t> alleles);

// get maximum supported ploidy
uint32_t get_max_genotype_ploidy();

// get maximum supported alleles
uint32_t get_max_genotype_alleles();


#endif // GENOTYPE_H
