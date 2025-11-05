# cython: language_level=3

from libcpp cimport bool
from libc.stdint cimport uint32_t, uint64_t
from . cimport cpp

cdef class Read:
	cdef cpp.Read *thisptr
	cdef bool ownsptr


cdef class ReadSet:
	cdef cpp.ReadSet *thisptr


cdef class GenotypeLikelihoods:
	cdef cpp.GenotypeLikelihoods *thisptr
	
	
cdef class Genotype:
	cdef cpp.Genotype *thisptr
	cdef uint64_t index


cdef class GenotypingAlgorithm:
	cdef cpp.GenotypingAlgorithm *thisptr