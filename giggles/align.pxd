# cython: language_level=3

from . cimport cpp

cdef class WFAWrapper:
	cdef cpp.WFAWrapper *thisptr