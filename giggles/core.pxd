# cython: language_level=3

from libcpp cimport bool
from libc.stdint cimport uint32_t, uint64_t
from . cimport cpp


cdef class NumericSampleIds:
	cdef dict mapping
	cdef bool frozen


cdef class Read:
	cdef cpp.Read *thisptr
	cdef bool ownsptr


cdef class ReadSet:
	cdef cpp.ReadSet *thisptr


cdef class Pedigree:
	cdef cpp.Pedigree *thisptr
	cdef NumericSampleIds numeric_sample_ids


cdef class PhredGenotypeLikelihoods:
	cdef cpp.PhredGenotypeLikelihoods *thisptr
	
	
cdef class Genotype:
	cdef cpp.Genotype *thisptr
	cdef uint64_t index
	cdef uint32_t ploidy


cdef class GenotypeHMM:
	cdef cpp.GenotypeHMM *thisptr
	cdef Pedigree pedigree
	cdef NumericSampleIds numeric_sample_ids
