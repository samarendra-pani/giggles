# cython: language_level=3
# Code modified from WhatsHap (https://github.com/whatshap/whatshap)
"""
Declarations for all external C++ classes that are wrapped from Cython.
"""

from libcpp.string cimport string
from libc.stdint cimport uint8_t

cdef extern from "../external/wrappers/wfawrapper.h":
	cdef cppclass WFAWrapper:
		WFAWrapper() except +
		int align(const char*, int, const char*, int, uint8_t) except +