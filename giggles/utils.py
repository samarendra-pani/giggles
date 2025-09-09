# Code taken from WhatsHap (https://github.com/whatshap/whatshap)

import gzip
import itertools
from typing import Sequence
import pyfaidx
from abc import ABC, abstractmethod

from giggles.logger import logger
from giggles import __version__
from giggles.core import (
    readselection,
    Genotype
)

def detect_file_format(path):
    """
    Detect file format and return 'BAM', 'CRAM', 'GAF' or None. None indicates an
    unrecognized file format.
    """
    with open(path, "rb") as f:
        first_bytes = f.read(16)
        if first_bytes.startswith(b"CRAM"):
            return "CRAM"

    gzip_header = b"\037\213"
    if first_bytes.startswith(gzip_header):
        with gzip.GzipFile(path, "rb") as f:
            first_bytes = f.read(16)
            if first_bytes.startswith(b"BAM\1"):
                return "BAM"
    
    #TODO: Need a better way to detect GAF file
    path = path.split(".")
    if "gaf" in path[-2:]:
        return "GAF"

    return None


def IndexedFasta(path):
    try:
        f = pyfaidx.Fasta(path, as_raw=True, sequence_always_upper=True, build_index=False)
    except pyfaidx.IndexNotFoundError:
        raise Exception(f"FASTA file {path} is not indexed")
    return f

def select_reads(readset, max_coverage, preferred_source_ids=None):
    logger.info(f"Reducing coverage to at most {max_coverage}X by selecting most informative reads ...")
    selected_indices = readselection(readset, max_coverage, preferred_source_ids)
    selected_reads = readset.subset(selected_indices)
    logger.info(f"Selected {len(selected_reads)} reads covering {len(selected_reads.get_positions())} variants")

    return selected_reads

class RecombinationCostComputer(ABC):
    @abstractmethod
    def compute(self, positions):
        pass

class UniformRecombinationCostComputer(RecombinationCostComputer):
    def __init__(self, recombination_rate, eff_pop_size):
        self._recombination_rate = recombination_rate
        self._eff_pop_size = eff_pop_size

    @staticmethod
    def uniform_recombination_map(recombrate, eff_pop_size, positions):

        # For a list of positions and a constant recombination rate (in cM/Mb),
        # return a list "results" of the same length as "positions" such that
        # results[i] is the phred-scaled recombination probability between
        # positions[i-1] and positions[i].
        
        return [(positions[i] - positions[i - 1])*recombrate*eff_pop_size*(4/(pow(10,6))) for i in range(1, len(positions))]

    def compute(self, positions):
        return self.uniform_recombination_map(self._recombination_rate, self._eff_pop_size, positions)


def bin_coeff(n, k):
    if (k < 0) or (n < 0) or (n < k):
        return 0
	
    result = 1.0
    if (k > n-k):
        k = n-k
	
    for i in range(k):
        result *= (n-i)
        result /= (i+1)
	
    return int(result)


def int_to_diploid_multiallelic_gt(numeric_repr):
    """Converts the classic numeric representation of multi-allelic, diploid genotypes
    into a genotype object
    """
    if numeric_repr == -1:
        return Genotype([])
    ploidy = 2
    genotype = [-1,-1]
    pth = ploidy
    max_allele_index = numeric_repr
    leftover_genotype_index = numeric_repr

    while (pth > 0):
        for allele_index in range(max_allele_index+1):
            i = bin_coeff(pth + allele_index - 1, pth)
            if (i >= leftover_genotype_index) or (allele_index == max_allele_index):
                if (i > leftover_genotype_index):
                    allele_index -= 1
                leftover_genotype_index -= bin_coeff(pth + allele_index - 1, pth)
                pth -= 1
                max_allele_index = allele_index
                genotype[pth] = allele_index
                break
    
    return Genotype(genotype)


def determine_genotype(likelihoods: Sequence[float], threshold_prob: float, n_allele: int) -> float:
    """given genotype likelihoods for 0/0, 0/1, 1/1, determines likeliest genotype"""

    assert bin_coeff(n_allele + 1, n_allele - 1) == len(likelihoods)
    to_sort = []
    for i in range(len(likelihoods)):
        to_sort.append((likelihoods[int_to_diploid_multiallelic_gt(i)], i))
    to_sort.sort(key=lambda x: x[0])

    # make sure there is a unique maximum which is greater than the threshold
    if (to_sort[-1][0] > to_sort[-2][0]) and (to_sort[-1][0]-to_sort[-2][0] > threshold_prob):
        return int_to_diploid_multiallelic_gt(to_sort[-1][1])
    else:
        return int_to_diploid_multiallelic_gt(-1)

def reverse_complement(seq):
    seq = seq.replace("A", "t").replace(
        "C", "g").replace("T", "a").replace("G", "c")
    seq = seq.upper()
    
    seq = seq[::-1]
    return seq

def reverse_cigar(cigar):
    all_cigars = ["".join(x) for _, x in itertools.groupby(cigar, key=str.isdigit)]
    new_cigar = ""
    for i in range(len(all_cigars), 0, -2):
        new_cigar += str(all_cigars[i - 2]) + str(all_cigars[i - 1])
    return new_cigar