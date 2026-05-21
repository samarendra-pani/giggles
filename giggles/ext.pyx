# cython: language_level=3

from libcpp.string cimport string
from libc.stdint cimport uint8_t, int32_t

"""
Wrapper for WFA2-lib
"""
cdef class WFAWrapper:
	def __cinit__(self, bandwidth):
		cdef int32_t _bandwidth = bandwidth
		self.thisptr = new extcpp.WFAWrapper(_bandwidth)
	
	def __dealloc__(self):
		del self.thisptr
	
	def align(self, text, pattern, state):
		cdef string _text = text.encode('UTF-8')
		cdef string _pattern = pattern.encode('UTF-8')
		cdef uint8_t _type = state

		return self.thisptr.align(_text, _pattern, _type)
