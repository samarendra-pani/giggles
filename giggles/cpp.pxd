# cython: language_level=3
# Code modified from WhatsHap (https://github.com/whatshap/whatshap)
"""
Declarations for all C++ classes that are wrapped from Cython.
"""
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libc.stdint cimport uint32_t, uint64_t
from libcpp.unordered_map cimport unordered_map


cdef extern from "../src/read.h":
	cdef cppclass Read:
		Read(string, int, int, int) except +
		Read(Read) except +
		vector[int] getMapqs() except +
		string getName() except +
		int getSourceID() except +
		int getReferenceStart() except +
		int getVariantCount() except +
		int getPosition(int) except +
		int getAllele(int) except +
		vector[unsigned int] getScores(int) except +
		void setPosition(int, int)  except +
		void setAllele(int, int) except +
		void setScores(int, vector[unsigned int]) except +
		void addVariant(int, int, vector[unsigned int]) except +
		void addHaplotag(string, int) except +
		void addMapq(int) except +
		void sortVariants() except +
		bool isSorted() except +
		

cdef extern from "../src/readset.h":
	cdef cppclass ReadSet:
		ReadSet() except +
		void add(Read*) except +
		string toString() except +
		int size() except +
		Read* get(int) except +
		Read* getByName(string, int) except +
		void sort() except +
		ReadSet* subset(IndexSet*) except +
		void assign_selection_status(IndexSet*) except +
		# TODO: Check why adding "except +" here doesn't compile
		vector[unsigned int]* get_positions() except +


cdef extern from "../src/genotypelikelihoods.h":
	cdef cppclass GenotypeLikelihoods:
		GenotypeLikelihoods(vector[long double], unsigned int) except +
		string toString() except +
		long double get_by_genotype(Genotype) except +
		void get_genotypes(vector[Genotype]&) except +
		unsigned int size() except +


cdef extern from "../src/binomial.h":
	cdef int binomial_coefficient(int n, int k) except +


cdef extern from "../src/genotype.h":
	cdef cppclass Genotype:
		Genotype() except +
		Genotype(vector[uint32_t]) except +
		vector[uint32_t] as_vector() except +
		bool is_none() except +
		uint64_t get_index() except +
		string toString() except +
	cdef bool operator==(Genotype,Genotype) except +
	cdef bool operator!=(Genotype,Genotype) except +
	cdef bool operator<(Genotype,Genotype) except +
	cdef vector[uint32_t] convert_index_to_alleles(uint64_t index) except +


cdef extern from "../src/indexset.h":
	cdef cppclass IndexSet:
		IndexSet() except +
		void add(int) except +


cdef extern from "../src/genotypingalgorithm.h":
	cdef cppclass GenotypingAlgorithm:
		GenotypingAlgorithm(ReadSet* readset, vector[float] recombcost, unsigned int n_samples, vector[unsigned int]* positions, vector[unsigned int]* n_allele_positions, vector[vector[int]]*, vector[bool]*) except +
		vector[long double] get_genotype_likelihoods(unsigned int position) except +