import gzip
import logging
from collections import defaultdict
from typing import Optional, DefaultDict
import pyfaidx

from giggles import __version__
from giggles.core import (
    readselection
)
from giggles.graph import ComponentFinder

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


def plural_s(n: int) -> str:
    return "" if n == 1 else "s"

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

def setup_families(samples):
    """
    Return families, family_trios pair.

    families maps a family representative to a list of family members

    family_trios maps a family representative to a list of trios in this family
    """

    # list of all trios across all families
    all_trios = dict()

    # Keep track of connected components (aka families) in the pedigree
    family_finder = ComponentFinder(samples)

    # map family representatives to lists of family members
    families = defaultdict(list)
    for sample in samples:
        families[family_finder.find(sample)].append(sample)

    # map family representatives to lists of trios for this family
    family_trios = defaultdict(list)
    for trio in all_trios:
        family_trios[family_finder.find(trio.child)].append(trio)
    logger.info(
        "Working on %d%s samples from %d famil%s",
        len(samples),
        plural_s(len(samples)),
        len(families),
        "y" if len(families) == 1 else "ies",
    )

    return families, family_trios
