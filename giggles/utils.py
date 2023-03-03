import gzip
import logging
from collections import defaultdict
from typing import Optional, DefaultDict

import pyfaidx
from dataclasses import dataclass


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
