import os
from typing import Optional, Iterable

import pysam
import heapq
from dataclasses import dataclass



@dataclass
class AlignmentWithSourceID:
    source_id: int
    bam_alignment: pysam.AlignedSegment


class AlignmentFileNotIndexedError(Exception):
    pass


class SampleNotFoundError(Exception):
    pass


class ReferenceNotFoundError(Exception):
    pass


class EmptyAlignmentFileError(Exception):
    pass

class SampleBamReader:
    """
    A wrapper for Samfile that provides only those reads from a BAM or CRAM file
    that belong to a specified sample. The BAM/CRAM file must have an index
    (bai/crai).
    """

    def __init__(self, path: str, *, source_id: int = 0, reference: Optional[str] = None):
        """
        path -- path to BAM/CRAM file
        reference -- path to FASTA reference
        """
        self.source_id: int = source_id
        if reference:
            reference = os.path.abspath(reference)

        self._samfile = pysam.AlignmentFile(path, reference_filename=reference)
        try:
            fetcher = self._samfile.fetch()
        except ValueError:
            raise AlignmentFileNotIndexedError(path)

        # For CRAM files, AlignmentFile.fetch() does not raise an error if the
        # reference could not be retrieved; it just returns an empty list.
        # To detect this situation, we treat all empty alignment files as an error.
        try:
            next(fetcher)
        except StopIteration:
            raise EmptyAlignmentFileError(path) from None
        self._references = frozenset(self._samfile.references)

    def has_reference(self, name: str) -> bool:
        return name in self._references

    def fetch(self, reference: str, start: int = 0, end: Optional[int] = None):
        """
        Yield instances of AlignmentWithSourceID, with source_id value given
        at construction time.
        Raise KeyError if sample not found among samples named in RG header.
        Raise ReferenceNotFoundError when reference is not in the BAM.
        """
        if reference not in self._references:
            # fetch() (below) raises a ValueError when the reference is invalid,
            # but also on other errors.
            raise ReferenceNotFoundError(reference)
        for bam_read in self._samfile.fetch(reference, start=start, stop=end):
            yield AlignmentWithSourceID(self.source_id, bam_read)
        
    def close(self) -> None:
        self._samfile.close()


class ComparableAlignedSegment:
    """
    Heapsort wants to be able to use the less than operator. Native
    AlignedSegment instances do not support this.
    """

    def __init__(self, aligned_segment: pysam.AlignedSegment, source_id: int):
        self.segment: pysam.AlignedSegment = aligned_segment
        self.source_id: int = source_id

    def __lt__(self, other) -> bool:
        self_pos = self.segment.reference_start
        other_pos = other.segment.reference_start
        return (self_pos < other_pos) or (
            self_pos == other_pos and self.source_id < other.source_id
        )


class MultiBamReader:
    """
    Read multiple sorted CRAM/BAM files and merge them on the fly.
    """

    def __init__(self, paths: Iterable[str], *, reference: Optional[str] = None):
        self._readers = []
        for source_id, path in enumerate(paths):
            self._readers.append(SampleBamReader(path, source_id=source_id, reference=reference))

    def fetch(
        self,
        reference: Optional[str] = None,
        start: int = 0,
        end: Optional[int] = None,
    ):
        """
        Yield reads from the specified region in all the opened CRAM/BAM files,
        merging them on the fly. Each CRAM/BAM file must be indexed.

        Yields instances of AlignmentWithSourceID, where source_id corresponds to
        index of BAM file name given at construction time.
        """
        assert reference is not None

        def make_comparable(reader):
            for alignment in reader.fetch(reference, start, end):
                yield ComparableAlignedSegment(alignment.bam_alignment, alignment.source_id)

        iterators = []
        for reader in self._readers:
            iterators.append(make_comparable(reader))
        for it in heapq.merge(*iterators):
            yield AlignmentWithSourceID(it.source_id, it.segment)

    def has_reference(self, name: str) -> bool:
        return all(reader.has_reference(name) for reader in self._readers)

    def close(self) -> None:
        for f in self._readers:
            f.close()