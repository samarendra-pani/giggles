# cython: language_level=3

from libcpp.string cimport string
from libc.stdint cimport uint8_t, int32_t

"""
Wrapper for WFA2-lib
"""
cdef class WFAWrapper:
	def __cinit__(self):
		self.thisptr = new extcpp.WFAWrapper()
	
	def __dealloc__(self):
		del self.thisptr
	
	def align(self, str text, str pattern, int state):
		# 1. Force explicit local bytes variables
		# Using typed variable definitions anchors their reference counts
		cdef bytes text_bytes = text.encode('UTF-8')
		cdef bytes pattern_bytes = pattern.encode('UTF-8')
		
		# 2. Safely cast to raw pointers
		cdef const char* c_text = text_bytes
		cdef int text_len = len(text_bytes)
		
		cdef const char* c_pattern = pattern_bytes
		cdef int pattern_len = len(pattern_bytes)
		
		cdef uint8_t _type = state

		# 3. Now it is 100% safe to execute. 
		# Python guarantees text_bytes and pattern_bytes live until this method returns.
		return self.thisptr.align(c_text, text_len, c_pattern, pattern_len, _type)
