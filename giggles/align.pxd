# cython: language_level=3

from . cimport extcpp

cdef class WFAWrapper:
	cdef extcpp.WFAWrapper *thisptr