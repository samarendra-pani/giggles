# cython: language_level=3

# Code modified from WhatsHap (https://github.com/whatshap/whatshap)

"""
Wrappers for core C++ classes.
"""
# Do not use the distutils directives here, but configure everything in
# setup.py, such as language and sources. It would work during development, but
# during a regular installation, the module will be compiled from the
# pre-generated .cpp file and the .pyx file is not read.

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport uint32_t, uint64_t
from . cimport cpp

from .variant import Variant
from collections import namedtuple
from cython.operator cimport dereference as deref

cdef class Read:
	def __cinit__(self, str name = None, int mapq = 0, int source_id = 0, int reference_start = -1):
		cdef string _name = b''
		if name is None:
			self.thisptr = NULL
			self.ownsptr = False
		else:
			# TODO: Is this the best way to handle string arguments?
			_name = name.encode('UTF-8')
			self.thisptr = new cpp.Read(_name, mapq, source_id, reference_start)
			self.ownsptr = True

	def __dealloc__(self):
		if self.ownsptr:
			assert self.thisptr != NULL
			del self.thisptr

	def __repr__(self):
		assert self.thisptr != NULL
		return 'Read(name={!r}, mapq={}, source_id={}, reference_start={}, variants={})'.format(
			self.name, self.mapqs, self.source_id, self.reference_start, list(self))

	property mapqs:
		def __get__(self):
			assert self.thisptr != NULL
			return tuple(self.thisptr.getMapqs())

	property name:
		def __get__(self):
			assert self.thisptr != NULL
			return self.thisptr.getName().decode('utf-8')

	property source_id:
		def __get__(self):
			assert self.thisptr != NULL
			return self.thisptr.getSourceID()
	
	property reference_start:
		def __get__(self):
			assert self.thisptr != NULL
			return self.thisptr.getReferenceStart()

	def __iter__(self):
		"""Iterate over all variants in this read"""
		assert self.thisptr != NULL
		for i in range(len(self)):
			yield self[i]

	def __len__(self):
		"""Return number of variants in this read"""
		assert self.thisptr != NULL
		return self.thisptr.getVariantCount()

	def __getitem__(self, key):
		"""Return Variant object at the given integer index"""
		assert self.thisptr != NULL
		if isinstance(key, slice):
			raise NotImplementedError("Read does not support slices")
		assert isinstance(key, int)
		cdef int n = self.thisptr.getVariantCount()
		if not (-n <= key < n):
			raise IndexError('Index out of bounds: {}'.format(key))
		if key < 0:
			key = n + key
		return Variant(
			position=self.thisptr.getPosition(key),
			allele=self.thisptr.getAllele(key),
			scores=self.thisptr.getScores(key),
		)

	def __setitem__(self, index, variant):
		assert self.thisptr != NULL
		cdef int n = self.thisptr.getVariantCount()
		if not (-n <= index < n):
			raise IndexError('Index out of bounds: {}'.format(index))
		if index < 0:
			index = n + index
		if not isinstance(variant, Variant):
			raise ValueError('Expected instance of Variant, but found {}'.format(type(variant)))
		self.thisptr.setPosition(index, variant.position)
		self.thisptr.setAllele(index, variant.allele)
		self.thisptr.setScores(index, variant.scores)

	def __contains__(self, position):
		"""Return whether this read contains a variant at the given position.
		A linear search is used.
		"""
		assert self.thisptr != NULL
		assert isinstance(position, int)
		for variant in self:
			if variant.position == position:
				return True
		return False
	
	def __getstate__(self):
		mapqs = [mapq for mapq in self.mapqs]
		variants = [(var.position, var.allele, var.scores) for var in self]
		return (mapqs, self.name, self.source_id, self.reference_start, variants)

	def __setstate__(self, state):
		mapqs, name, source_id, reference_start, variants = state

		# TODO: Duplicated code from __cinit__ is ugly, but cinit cannot be used here directly
		cdef string _name = b''
		if name is None:
			self.thisptr = NULL
			self.ownsptr = False
		else:
			# TODO: Is this the best way to handle string arguments?
			_name = name.encode('UTF-8')
			self.thisptr = new cpp.Read(_name, mapqs[0] if len(mapqs) > 0 else 0, source_id, reference_start)
			self.ownsptr = True

		for mapq in mapqs[1:]:
			self.add_mapq(mapq)
		for (pos, allele, scores) in variants:
			self.add_variant(pos, allele, scores)

	def add_variant(self, int position, int allele, vector[unsigned int] scores):
		assert self.thisptr != NULL
		cdef vector[unsigned int] c_scores
		c_scores.resize(len(scores))
		for i in range(len(scores)):
			c_scores[i] = scores[i]
		self.thisptr.addVariant(position, allele, c_scores)

	def add_haplotag(self, str hp, int ps):
		cdef string _hp = b''
		_hp = hp.encode('UTF-8')
		self.thisptr.addHaplotag(_hp, ps)

	def add_mapq(self, int mapq):
		assert self.thisptr != NULL
		self.thisptr.addMapq(mapq)

	def sort(self):
		assert self.thisptr != NULL
		self.thisptr.sortVariants()

	def is_sorted(self):
		assert self.thisptr != NULL
		return self.thisptr.isSorted()


cdef class ReadSet:
	def __cinit__(self):
		self.thisptr = new cpp.ReadSet()

	def __dealloc__(self):
		del self.thisptr

	def add(self, Read read):
		"""Adds a read to the set.
		WARNING: this will internally create a copy of the wrapped C++ Read object,
		so that subsequent changes to the Read don't affect the
		newly created copy that is added to the ReadSet."""
		self.thisptr.add(new cpp.Read(read.thisptr[0]))

	def __str__(self):
		return self.thisptr.toString().decode('utf-8')

	def __iter__(self):
		for i in range(self.thisptr.size()):
			yield self[i]

	def __len__(self):
		return self.thisptr.size()

	def __getitem__(self, key):
		if isinstance(key, slice):
			raise NotImplementedError('ReadSet does not support slices')
		cdef string name = b''
		cdef cpp.Read* cread = NULL
		cdef Read read = Read()
		if isinstance(key, int):
			read.thisptr = self.thisptr.get(key)
		elif isinstance(key, str):
			raise NotImplementedError('Querying a ReadSet by read name is deprecated, please query by (source_id, name) instead')
		elif isinstance(key, tuple) and (len(key) == 2) and (isinstance(key[0],int) and isinstance(key[1],str)):
			source_id = key[0]
			name = key[1].encode('UTF-8')
			cread = self.thisptr.getByName(name, source_id)
			if cread == NULL:
				raise KeyError(key)
			else:
				read.thisptr = cread
		else:
			assert False, 'Invalid key: {}'.format(key)
		return read
	
	def __getstate__(self):
		return ([read for read in self])
	
	def __setstate__(self, state):
		self.thisptr = new cpp.ReadSet()
		for read in state:
			self.add(read)

	def sort(self):
		"""Sort contained reads by the position of the first variant they contain. Note that
		this is not necessarily the variant with the lowest position, unless sort() has been
		called on all contained reads. Ties are resolved by comparing the read name."""
		self.thisptr.sort()

	def subset(self, reads_to_select):
		# TODO: is there a way of avoiding to unecessarily creating/destroying a ReadSet object?
		cdef cpp.IndexSet* index_set = new cpp.IndexSet()
		cdef int i
		for i in reads_to_select:
			index_set.add(i)
		result = ReadSet()
		del result.thisptr
		result.thisptr = self.thisptr.subset(index_set)
		del index_set
		return result
	
	def assign_selection_status(self, reads_to_select):
		cdef cpp.IndexSet* index_set = new cpp.IndexSet()
		cdef int i
		for i in reads_to_select:
			index_set.add(i)
		self.thisptr.assign_selection_status(index_set)
		del index_set

	def get_positions(self):
		cdef vector[unsigned int]* v = self.thisptr.get_positions()
		result = list(v[0])
		del v
		return result

cdef class GenotypeLikelihoods:
	def __cinit__(self, vector[long double] gl, unsigned int num_alleles):
		self.thisptr = new cpp.GenotypeLikelihoods(gl, num_alleles)

	def __dealloc__(self):
		del self.thisptr

	def __str__(self):
		return self.thisptr.toString().decode('utf-8')

	def __getitem__(self, Genotype genotype):
		assert self.thisptr != NULL
		return self.thisptr.get_by_genotype(genotype.thisptr[0])

	def __len__(self):
		return self.thisptr.size()

	def __iter__(self):
		for genotype in self.genotypes():
			yield self[genotype]
			
	def __eq__(self, GenotypeLikelihoods other):
		if self.genotypes() != other.genotypes():
			return False
		for genotype in self.genotypes():
			if self[genotype] != other[genotype]:
				return False
		return True
		
	def genotypes(self):
		cdef vector[cpp.Genotype]* genotypes = new vector[cpp.Genotype]()
		self.thisptr.get_genotypes(deref(genotypes))
		result = [Genotype(genotype.as_vector()) for genotype in genotypes[0]]
		del genotypes
		return result
	
	
def binomial_coefficient(int n, int k):
	return cpp.binomial_coefficient(n, k)

			
cdef class Genotype:
	def __cinit__(self, vector[uint32_t] alleles):
		self.thisptr = new cpp.Genotype(alleles)
		self.index = self.thisptr.get_index()

	def __dealloc__(self):
		del self.thisptr

	def __str__(self):
		return self.thisptr.toString().decode('utf-8')
	
	def __repr__(self):
		return self.thisptr.toString().decode('utf-8')

	def is_none(self):
		return self.thisptr.is_none()

	def get_index(self):
		return self.thisptr.get_index()

	def as_vector(self):
		result = []
		cdef vector[uint32_t] alleles = self.thisptr.as_vector()
		for allele in alleles:
			result.append(allele)
		return alleles

	def __eq__(self, Genotype g):
		return self.thisptr[0] == g.thisptr[0]

	def __ne__(self, Genotype g):
		return self.thisptr[0] != g.thisptr[0]
	
	def __lt__(self, Genotype g):
		return self.thisptr[0] < g.thisptr[0]
	
	def __hash__(self):
		return hash(tuple(self.alleles))
	
	def __reduce__(self):
		# a tuple as specified in the pickle docs - (class_or_constructor, 
		# (tuple, of, args, to, constructor))
		cdef vector[uint32_t] alleles = cpp.convert_index_to_alleles(self.index)
		return (self.__class__, tuple([alleles]))
	


cdef class GenotypingAlgorithm:
	def __cinit__(self, ReadSet readset, recombcost, n_haplotypes, positions, n_allele_positions, allele_references, is_sv):
		"""
		The GenotypingAlgorithm performs an iterative phasing-genotyping algorithm
		using the following ReadSet at positions specified.
		"""
		# Prepare C++ vectors for optional arguments
		cdef vector[unsigned int] c_positions_stack
		cdef vector[unsigned int] c_n_allele_positions_stack
		cdef vector[vector[int]] c_allele_references_stack
		cdef vector[bool] c_is_sv_stack
		
		# Prepare pointers for optional arguments
		cdef vector[unsigned int]* c_positions_ptr = NULL
		cdef vector[unsigned int]* c_n_allele_positions_ptr = NULL
		cdef vector[vector[int]]* c_allele_references_ptr = NULL
		cdef vector[bool]* c_is_sv_ptr = NULL
		
		cdef unsigned int n_references = n_haplotypes
		
		# Fill C++ vectors and set pointers
		if positions is not None:
			for pos in positions:
				c_positions_stack.push_back(pos)
			c_positions_ptr = &c_positions_stack
		if n_allele_positions is not None:
			for pos in n_allele_positions:
				c_n_allele_positions_stack.push_back(pos)
			c_n_allele_positions_ptr = &c_n_allele_positions_stack
		if allele_references is not None:
			c_allele_references_stack.resize(len(allele_references))
			for ix, pos in enumerate(allele_references):
				for hap in pos:
					c_allele_references_stack.at(ix).push_back(hap)
			c_allele_references_ptr = &c_allele_references_stack
		if is_sv is not None:
			for sv in is_sv:
				c_is_sv_stack.push_back(sv)
			c_is_sv_ptr = &c_is_sv_stack
		
		# Finally, create the C++ object
		self.thisptr = new cpp.GenotypingAlgorithm(readset.thisptr, recombcost, n_references, c_positions_ptr, c_n_allele_positions_ptr, c_allele_references_ptr, c_is_sv_ptr)

	def __dealloc__(self):
		del self.thisptr

	def get_genotype_likelihoods(self, unsigned int pos, unsigned int num_allele):
		return GenotypeLikelihoods(self.thisptr.get_genotype_likelihoods(pos), num_alleles = num_allele)


include 'readselect.pyx'