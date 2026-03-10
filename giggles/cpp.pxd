# cython: language_level=3
# Code modified from WhatsHap (https://github.com/whatshap/whatshap)
"""
Declarations for all C++ classes that are wrapped from Cython.
"""
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libc.stdint cimport uint32_t, uint8_t
from libcpp.unordered_map cimport unordered_map


cdef extern from "../src/read.h":
	cdef cppclass Read:
		Read(string, uint32_t, uint32_t) except +
		Read(Read) except +
		vector[uint32_t] getMapqs() except +
		string getName() except +
		uint32_t getSourceID() except +
		uint32_t getVariantCount() except +
		uint32_t getPosition(uint32_t) except +
		vector[long double] getEmissionScores(uint32_t) except +
		void setPosition(uint32_t, uint32_t)  except +
		void setScores(uint32_t, vector[uint32_t]) except +
		void setEmissionScores(uint32_t, vector[long double]) except +
		void addVariant(uint32_t, vector[uint32_t]) except +
		void addVariant(uint32_t, vector[long double]) except +
		void setHaplotag(string) except +
		void setPhaseSet(int) except +
		void addMapq(uint32_t) except +
		void sortVariants() except +
		bool isSorted() except +
		

cdef extern from "../src/readset.h":
	cdef cppclass ReadSet:
		ReadSet() except +
		void add(Read*) except +
		string toString() except +
		uint32_t size() except +
		Read* get(uint32_t) except +
		Read* getByName(string, int) except +
		void sort() except +
		ReadSet* subset(IndexSet*) except +
		void assign_selection_status(IndexSet*) except +
		# TODO: Check why adding "except +" here doesn't compile
		vector[uint32_t]* get_positions() except +


cdef extern from "../src/genotypelikelihoods.h":
	cdef cppclass GenotypeLikelihoods:
		GenotypeLikelihoods(vector[long double], uint32_t, uint32_t) except +
		string toString() except +
		long double get_by_index(uint32_t) except +
		uint32_t size() except +


cdef extern from "../src/binomial.h":
	cdef int binomial_coefficient(int n, int k) except +


cdef extern from "../src/genotype.h":
	cdef cppclass Genotype:
		Genotype() except +
		Genotype(vector[uint32_t]) except +
		vector[uint32_t] as_vector() except +
		uint32_t get_ploidy() except +
		bool is_none() except +
		uint32_t get_index() except +
		string toString() except +
	cdef bool operator==(Genotype,Genotype) except +
	cdef bool operator!=(Genotype,Genotype) except +
	cdef bool operator<(Genotype,Genotype) except +
	cdef vector[uint32_t] convert_index_to_alleles(uint32_t index, uint32_t ploidy) except +
	cdef uint32_t get_max_genotype_ploidy() except +
	cdef uint32_t get_max_genotype_alleles() except +


cdef extern from "../src/indexset.h":
	cdef cppclass IndexSet:
		IndexSet() except +
		void add(int) except +


cdef extern from "../src/genotypingalgorithm.h":
	cdef cppclass GenotypingAlgorithm:
		GenotypingAlgorithm(ReadSet* readset, vector[float] recombcost, uint32_t n_samples, uint32_t ploidy, vector[uint32_t]* positions, vector[uint32_t]* n_allele_positions, vector[vector[int]]*, vector[bool]*) except +
		vector[long double] get_genotype_likelihoods(uint32_t position) except +

cdef extern from "../src/tests.h":
	cdef void run_all_tests() except +