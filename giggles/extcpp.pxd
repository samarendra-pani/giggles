# cython: language_level=3
# Code modified from WhatsHap (https://github.com/whatshap/whatshap)
"""
Declarations for all external C++ classes that are wrapped from Cython.
"""

from libcpp.string cimport string
from libc.stdint cimport uint32_t, uint8_t


cdef extern from "../external/wrappers/wfawrapper.h":
	cdef cppclass WFAWrapper:
		WFAWrapper(uint32_t) except +
		uint32_t align(string, string, uint8_t) except +