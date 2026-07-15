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
	
	def align(self, str allele, str query, int state):
		cdef bytes allele_bytes = allele.encode('UTF-8')
		cdef bytes query_bytes = query.encode('UTF-8')
		
		cdef const char* c_allele = allele_bytes
		cdef int allele_len = len(allele_bytes)
		
		cdef const char* c_query = query_bytes
		cdef int query_len = len(query_bytes)
		
		cdef uint8_t _type = state

		if state == 0:
			if allele_len >= 1.2*query_len or allele_len <= query_len/1.2:
				raise RuntimeError('Allele and Query lengths mismatch.')
		elif state != 3:
			if allele_len <= query_len/1.2:
				raise RuntimeError('Allele and Query lengths mismatch.')

		return self.thisptr.align(c_allele, allele_len, c_query, query_len, _type)
