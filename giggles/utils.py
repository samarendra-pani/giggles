# Code taken from WhatsHap (https://github.com/whatshap/whatshap)

import gzip
import logging
from collections import defaultdict
from typing import DefaultDict
import pyfaidx
from abc import ABC, abstractmethod

from giggles import __version__
from giggles.core import (
    readselection
)

class FastaNotIndexedError(Exception):
    pass


class InvalidRegion(Exception):
    pass


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
        raise FastaNotIndexedError(path)
    return f


_warning_count: DefaultDict[str, int] = defaultdict(int)


def warn_once(logger, msg: str, *args) -> None:
    if _warning_count[msg] == 0 and not logger.isEnabledFor(logging.DEBUG):
        logger.warning(msg + " Hiding further warnings of this type, use --debug to show", *args)
    else:
        logger.debug(msg, *args)
    _warning_count[msg] += 1


logger = logging.getLogger(__name__)


def select_reads(readset, max_coverage, preferred_source_ids=None):
    logger.info(
        "Reducing coverage to at most %dX by selecting most informative reads ...", max_coverage
    )
    selected_indices = readselection(readset, max_coverage, preferred_source_ids)
    selected_reads = readset.subset(selected_indices)
    logger.info(
        "Selected %d reads covering %d variants",
        len(selected_reads),
        len(selected_reads.get_positions()),
    )

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