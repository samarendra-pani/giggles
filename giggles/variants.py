"""
Detect variants in reads.
"""
import logging
import math
from collections import defaultdict, Counter, namedtuple
from typing import Iterable, Iterator, List, Optional
import re
from pywfa import WavefrontAligner

from giggles.core import Read, ReadSet, NumericSampleIds
from giggles.bam import SampleBamReader, MultiBamReader, BamReader
from giggles.gaf import GafParser, SampleGafParser
from giggles.align import edit_distance, edit_distance_affine_gap
from giggles._variants import _iterate_cigar

logger = logging.getLogger(__name__)

class CommandLineError(Exception):
    """An anticipated command-line error occurred. This ends up as a user-visible error message"""

class ReadSetError(Exception):
    pass


class AlignmentReader:
    """
    Superclass for the GAF and BAM Readers
    """
    def __init__(
            self,
            paths,
            numeric_sample_ids: NumericSampleIds,
            mapq_threshold: int,
            realign_mode: str,
            overhang: int,
            gap_start: int,
            gap_extend: int,
            default_mismatch: int,
            em_prob_params: List[float]):

        self._paths = paths
        self._mapq_threshold = mapq_threshold
        self._numeric_sample_ids = numeric_sample_ids
        self._realign_mode = realign_mode
        if realign_mode == "ed":
            self._aligner = edit_distance
        elif realign_mode == "wfa_full":
            self._aligner = WavefrontAligner(mismatch=default_mismatch, 
                                         gap_opening=gap_start,
                                         gap_extension=gap_extend)
        elif realign_mode == "wfa_score":
            self._aligner = WavefrontAligner(mismatch=default_mismatch, 
                                         gap_opening=gap_start,
                                         gap_extension=gap_extend,
                                         scope='score')
        self._gap_start = gap_start
        self._gap_extend = gap_extend
        self._default_mismatch = default_mismatch
        self._overhang = overhang
        self._em_params = em_prob_params
        
    @property
    def n_paths(self):
        return len(self._paths)

    @staticmethod
    def _make_readset_from_grouped_reads(groups: Iterable[List[Read]]) -> ReadSet:
        read_set = ReadSet()
        for group in groups:
            read_set.add(merge_reads(*group))
        return read_set

    @staticmethod
    def split_cigar(cigar, i, consumed):
        """
        Split a CIGAR into two parts. i and consumed describe the split position.
        i is the element of the cigar list that should be split, and consumed says
        at how many operations to split within that element.

        The CIGAR is given as a list of (operation, length) pairs.

        i -- split at this index in cigar list
        consumed -- how many cigar ops at cigar[i] are to the *left* of the
            split position

        Return a tuple (left, right).

        Example:
        Assume the cigar is 3M 1D 6M 2I 4M.
        With i == 2 and consumed == 5, the cigar is split into
        3M 1D 5M and 1M 2I 4M.
        """
        middle_op, middle_length = cigar[i]
        assert consumed <= middle_length
        if consumed > 0:
            left = cigar[:i] + [(middle_op, consumed)]
        else:
            left = cigar[:i]
        if consumed < middle_length:
            right = [(middle_op, middle_length - consumed)] + cigar[i + 1 :]
        else:
            right = cigar[i + 1 :]
        return left, right

    @staticmethod
    def cigar_prefix_length(cigar, reference_bases):
        """
        Given a prefix of length reference_bases relative to the reference, how
        long is the prefix of the read? In other words: If reference_bases on
        the reference are consumed, how many bases on the query does that
        correspond to?

        If the position is within or at the end of an insertion (which do not
        consume bases on the reference), then the number of bases up to the
        beginning of the insertion is reported.

        Return a pair (reference_bases, query_bases) where the value for
        reference_bases may be smaller than the requested one if the CIGAR does
        not cover enough reference bases.

        Reference skips (N operators) are treated as the end of the read. That
        is, no positions beyond a reference skip are reported.
        """
        ref_pos = 0
        query_pos = 0
        for op, length in cigar:
            if op in (0, 7, 8):  # M, X, =
                ref_pos += length
                query_pos += length
                if ref_pos >= reference_bases:
                    return (reference_bases, query_pos + reference_bases - ref_pos)
            elif op == 2:  # D
                ref_pos += length
                if ref_pos >= reference_bases:
                    return (reference_bases, query_pos)
            elif op == 1:  # I
                query_pos += length
            elif op == 4 or op == 5:  # soft or hard clipping
                pass
            elif op == 3:  # N
                # Always stop at reference skips
                return (reference_bases, query_pos)
            else:
                assert False, "unknown CIGAR operator"
        assert ref_pos < reference_bases
        return (ref_pos, query_pos)

    @staticmethod
    def realign(
            aligner: WavefrontAligner,
            variant,
            bam_read,
            cigartuples,
            i,
            consumed,
            query_pos,
            reference,
            mode,
            overhang,
            emission_parameters):
        """
        Realign a read to the two alleles of a single variant.
        i and consumed describe where to split the cigar into a part before the
        variant position and into a part starting at the variant position, see split_cigar().

        variant -- VcfVariant
        bam_read -- the AlignedSegment
        cigartuples -- the AlignedSegment.cigartuples property (accessing it is expensive, so re-use it)
        i, consumed -- see split_cigar method
        query_pos -- index of the query base that is at the variant position
        reference -- the reference as a str-like object (unlike original implementation, this is only the sequence of the alignment path and not the whole chromosome)
        overhang -- extend alignment by this many bases to left and right
        gap_start, gap_extend -- use these parameters for affine gap cost alignment
        default_mismatch -- use this as mismatch cost in case no base qualities are in alignment
        emission_parameters -- a list which contains the probabilities of the cigar being match, mismatch, insertion, and deletion.
        """
        # Do not process symbolic alleles like <DEL>, <DUP>, etc.
        if any([alt.startswith("<") for alt in variant.alternative_allele]):
            return None, None

        # There is a big difference between the previous implementation and what is needed.
        # In the previous code, the CIGAR is against the reference always and hence we need to realign only for the alternate alleles.
        # With GAF, the CIGAR is not always against the reference (sometimes it is not against ref or any of the alt and can be with a path that is not an allele traversal)
        # So we need to generalize the process to realign using the variant record and the cigar tuples.

        left_cigar, right_cigar = AlignmentReader.split_cigar(cigartuples, i, consumed)

        left_ref_bases, left_query_bases = AlignmentReader.cigar_prefix_length(
            left_cigar[::-1], overhang
        )
        
        if variant.reference_allele == "*":
            ref_allele = ""
        else:
            ref_allele = variant.reference_allele

        right_ref_bases, right_query_bases = AlignmentReader.cigar_prefix_length(
            right_cigar, len(ref_allele) + overhang
        )

        assert variant.position - left_ref_bases >= 0
        assert variant.position + right_ref_bases <= len(reference)

        query = bam_read.query_sequence[
            query_pos - left_query_bases : query_pos + right_query_bases
        ]
        
        left_overhang = reference[variant.position - left_ref_bases : variant.position]
        right_overhang = reference[variant.position + right_ref_bases - overhang : variant.position + right_ref_bases]
        
        ref = left_overhang + ref_allele + right_overhang
        
        alts = []
        for alt_allele in variant.alternative_allele:
            if alt_allele != "*":
                alt = left_overhang + alt_allele + right_overhang
            else:
                alt = left_overhang + right_overhang
            alts.append(alt)
        
        prob = []
        min_prob = 0
        min_allele = None
        for index, allele in enumerate([ref]+alts):
            if mode == "ed":
                prob.append(-aligner(query, allele))    #edit distance is positive. Need to change to negative.
            elif mode == "wfa_score":
                prob.append(aligner(query, allele).score)
            elif mode == "wfa_full":
                prob.append(AlignmentReader.calculate_emission_log_probability(aligner(query, allele).cigartuples, emission_parameters))
            if prob[index] < min_prob:
                min_prob = prob[index]
                min_allele = index
        
        base_qual_score = 30
        
        return min_allele, prob, base_qual_score
        
    @staticmethod
    def calculate_emission_log_probability(cg, params):
        """
        CIGAR Operation to Number Conversion (https://github.com/kcleal/pywfa/tree/master):
        M: 0
        I: 1
        D: 2
        N: 3
        S: 4
        H: 5
        =: 7
        X: 8
        B: 9

        The output cigar considers M as match.
        """
        prob = 0
        count = {'M': 0, 'X': 0, 'I': 0, 'D': 0}
        for op,c in cg:
            if op == 0:
                count['M'] += c
            elif op == 1:
                count['I'] += c
            elif op == 2:
                count['D'] += c
            elif op == 8:
                count['X'] += c
        for i, op in enumerate(count.keys()):
            prob += count[op]*math.log(params[i])
        return prob

    @staticmethod
    def detect_alleles_by_alignment(
        aligner,
        variants,
        j,
        read,
        reference,
        mode,
        overhang=10,
        emission_parameters=None
    ):
        """
        Detect which alleles the given bam_read covers. Detect the correct
        alleles of the variants that are covered by the given bam_read.

        Yield tuples (position, allele, quality).

        variants -- list of variants (VcfVariant objects)
        j -- index of the first variant (in the variants list) to check
        """
        # Accessing bam_read.cigartuples is expensive, do it only once
        cigartuples = read.cigartuples

        # For the same reason, the following check is here instad of
        # in the _usable_alignments method
        if not cigartuples:
            return
        for index, i, consumed, query_pos in _iterate_cigar(variants, j, read, cigartuples):
            allele, emission, quality = AlignmentReader.realign(
                aligner,
                variants[index],
                read,
                cigartuples,
                i,
                consumed,
                query_pos,
                reference,
                mode,
                overhang,
                emission_parameters
            )
            
            if allele is not None:
                yield (index, allele, emission, quality)


class GAFReader(AlignmentReader):
    """
    Associate VCF variant with GAF Read.
    """

    def __init__(
        self,
        paths: List[str],
        reference: str,
        read_fasta: str,
        numeric_sample_ids: NumericSampleIds,
        mapq_threshold: int = 20,
        realign_mode: str = "ed",
        overhang: int = 10,
        gap_start: int = 3,
        gap_extend: int = 1,
        default_mismatch: int = 2,
        em_prob_params: List[float] = [0.85, 0.05, 0.05, 0.05]
    ):
        super().__init__(paths, numeric_sample_ids, mapq_threshold, realign_mode, overhang, gap_start, gap_extend, default_mismatch, em_prob_params)
        self._reader: GafParser
        if len(paths) == 1:
            self._reader = SampleGafParser(paths[0], reference=reference, read_fasta=read_fasta, mapq=self._mapq_threshold)
        else:
            raise CommandLineError("Giggles does not support multiple GAF file parsing. Please provide a single GAF file.")

    def has_reference(self, chromosome):
        return self._reader.has_reference(chromosome)

    def read(self, chromosome, variants, sample=None, reference=None) -> ReadSet:
        """
        Detect alleles and return a ReadSet object containing reads representing
        the given variants.

        Using the provided reference, re-alignment is done which generates scores
        to be used in the HMM.

        chromosome -- name of chromosome to work on
        variants -- list of vcf.VcfVariant objects
        sample -- name of sample to work on. If None, read group information is
            ignored and all reads in the file are used.
        reference -- reference sequence of the given chromosome (or None)
        regions -- list of start,end tuples (end can be None)
        """
        # Since variants are identified by position, positions must be unique.
        if __debug__ and variants:
            varposc = Counter(variant.position for variant in variants)
            pos, count = varposc.most_common()[0]
            assert count == 1, "Position {} occurs more than once in variant list.".format(pos)

        logger.debug("Extracting Usable Alignments")
        alignments = self._usable_alignments(chromosome, variants)
        logger.debug("Converting Alignments to Read Objects")
        reads = self._alignments_to_reads(alignments, variants, sample)
        logger.debug("Grouping Reads into ReadSet Object")
        grouped_reads = self._remove_duplicate_reads(reads)
        logger.debug("ReadSet Object Successfully Created")
        readset = self._make_readset_from_grouped_reads(grouped_reads)
        return readset      

    @staticmethod
    def _remove_duplicate_reads(reads: Iterable[Read]) -> Iterator[List[Read]]:
        """
        remove reads which have been mapped multiple times and select one best read
        """
        groups = defaultdict(list)
        for read in reads:
            if groups[(read.source_id, read.name, read.sample_id)] == []:
                groups[(read.source_id, read.name, read.sample_id)] = [read]       # Keeping this as a list so that I dont need to change _make_readset_from_grouped_reads()
            else:
                old_read = groups[(read.source_id, read.name, read.sample_id)][0]

                # Check the number of variants it covers
                if len(old_read) < len(read):
                    groups[(read.source_id, read.name, read.sample_id)] = [read]
                
                # Check the mapping quality
                if old_read.mapqs < read.mapqs:
                    groups[(read.source_id, read.name, read.sample_id)] = [read]

        for group in groups.values():
            if len(group) > 1:
                raise ReadSetError(
                    "Read name {!r} occurs more than twice in the input file".format(group[0].name)
                )
            yield group
    
    @staticmethod
    def reverse_complement(seq):
        seq = seq.replace("A", "t").replace(
            "C", "g").replace("T", "a").replace("G", "c")
        seq = seq.upper()
        
        seq = seq[::-1]
        return seq


    def _usable_alignments(self, chromosome, variants):
        """
        Retrieve usable (suficient mapping quality, not secondary etc.)
        alignments from the alignment file
        """
        
        for alignment in self._reader(contig=chromosome):
            yield alignment

    def _alignments_to_reads(self, alignments, variants, sample):
        """
        Convert GAF alignments to Read objects.

        If reference is not None, alleles are detected through re-alignment.

        Yield Read objects.
        """
        numeric_sample_id = 0 if sample is None else self._numeric_sample_ids[sample]
        rgfa = self._reader._reference
        id_to_index = {}
        for i, v in enumerate(variants):
            id = v.id.split('>')
            start = id[1]
            end = id[-1]
            id_to_index[start] = [end, i] # Creating a look up table to find bubble end and index of variant from bubble start
        
        cg_letter_to_op = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, 'X': 7, '=': 8}
        Alignment = namedtuple('Alignment', ['cigartuples', 'reference_start', 'query_sequence'])        # Class created to maintain compatibility with old code
        for alignment in alignments:
            # Checking if the alignment is in reverse direction.
            
            reverse = False
            orient  = None
            for n in alignment.path:
                if n in ['>', '<']:
                    orient = n
                    continue
                if rgfa.get_node(n).tags["NO"] == 0:
                    reverse = orient == '<'
                    break
            
            # Reversing path and updating other variables if required
            if reverse:
                reverse_orient = {'>': '<', '<': '>'}
                new_alignment = []
                for n in alignment.path[::-1]:
                    if n in ['>', '<']:
                        new_alignment.insert(-1, reverse_orient[n])
                    else:
                        new_alignment.append(n)
                qs = alignment.q_start
                qe = alignment.q_end
                ql = alignment.q_len
                ps = alignment.p_start
                pe = alignment.p_end
                pl = alignment.p_len
                
                alignment.q_start = ql - qe
                alignment.q_end = ql - qs
                alignment.orient = '+'
                alignment.p_start = pl - pe
                alignment.p_end = pl - ps
                alignment.sequence = GAFReader.reverse_complement(alignment.sequence)
                alignment.path = new_alignment

            # Process the alignment path and change variant positions (which are now based on the paths).
            # The new definition of variant position is only useful for CIGAR processing. For all the downstream processing (especially for sorting), we need the actual position in the reference backbone.
            # Also output the reference sequence (which the sequence of the path)
            reference = ""
            variants_in_alignment = []
            end_found = False
            end = None
            for n in alignment.path:
                if n in ['>', '<']:
                    orient = n
                    continue
                node = rgfa.get_node(n)
                # Reference updated
                if orient == '<':
                    reference += GAFReader.reverse_complement(node.sequence)
                else:
                    reference += node.sequence
                # Finding variants and variant positions
                
                if node.tags['NO'] == 0:
                    try:
                        index = id_to_index[n][1]
                        if end != None:
                            end_found = end == n
                        end = id_to_index[n][0]
                    except KeyError:
                        continue
                    v = variants[index]
                    v.position = len(reference)
                    variants_in_alignment.append(v)
            
            if (len(variants_in_alignment) == 0) or (len(variants_in_alignment) == 1 and not end_found):
                continue

            if not end_found:
                variants_in_alignment.pop(-1)


            # Extract the aligned segement from the complete read sequence and create a new object.
            # Need cigartuples, and reference_start (where it starts in the reference. So the path start in this case.)


            gaf_aligned_segment = alignment.sequence[alignment.q_start:alignment.q_end]
            cg_tuples = []
            cg = list(filter(None, re.split("([MIDNSHP=X])", alignment.cigar)))
            for i in range(0,len(cg),2):
                l = int(cg[i])
                op = cg_letter_to_op[cg[i+1]]
                cg_tuples.append((op,l))
            
            # This new variable is created to make the gaf alignments compatible with the old code.
            processed_alignment = Alignment(cigartuples=cg_tuples, reference_start=alignment.p_start, query_sequence=gaf_aligned_segment)
            
            start_on_ref = rgfa.get_node(alignment.path[1]).start + alignment.p_start

            barcode = ""
            if alignment.has_tag("BX"):
                barcode = alignment.get_tag("BX")

            read = Read(
                alignment.read_id,
                alignment.mapping_quality,
                alignment.source_id,
                numeric_sample_id,
                start_on_ref,
                barcode,
            )

            detected = self.detect_alleles_by_alignment(
                self._aligner,
                variants_in_alignment,
                0,                              # Here this has been hardcoded to 0. In original code, this was the index of the first variant (index in the big list of variants) in the read. But now we have new list of variants just for this alignment.
                processed_alignment,
                reference,
                self._realign_mode,
                self._overhang,
                self._em_params)
        
            for j, allele, em, quality in detected:
                read.add_variant(variants_in_alignment[j].position_on_ref, allele, em, quality)
            if read:  # At least one variant covered and detected
                yield read

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        self._aligner.__dealloc__()
        self._reader.close()

class ReadSetReader(AlignmentReader):
    """
    Associate VCF variants with BAM reads.

    A VCF file contains variants, and a BAM file contain reads, but the
    information which read contains which variant is not available. This
    class re-discovers the variants in each read, using the
    knowledge in the VCF of where they should occur.
    """

    def __init__(
        self,
        paths: List[str],
        reference: Optional[str],
        numeric_sample_ids: NumericSampleIds,
        mapq_threshold: int = 20,
        realign_mode: str = "ed",
        overhang: int = 10,
        gap_start: int = 3,
        gap_extend: int = 1,
        default_mismatch: int = 2,
        em_prob_params: List[float] = [0.85, 0.05, 0.05, 0.05]
    ):
        """
        paths -- list of BAM paths
        reference -- path to reference FASTA (can be None)
        numeric_sample_ids -- sample ids in numeric format
        mapq_threshold -- minimum mapping quality
        overhang -- extend alignment by this many bases to left and right
        gap_start, gap_extend, default_mismatch -- parameters for affine gap cost alignment
        """
        super().__init__(paths, numeric_sample_ids, mapq_threshold, realign_mode, overhang, gap_start, gap_extend, default_mismatch, em_prob_params)
        self._reader: BamReader
        if len(paths) == 1:
            self._reader = SampleBamReader(paths[0], reference=reference)
        else:
            self._reader = MultiBamReader(paths, reference=reference)

    def has_reference(self, chromosome):
        return self._reader.has_reference(chromosome)

    def read(self, chromosome, variants, sample, reference) -> ReadSet:
        """
        Detect alleles and return a ReadSet object containing reads representing
        the given variants.

        If a reference is provided (reference is not None), alleles are
        detected by re-aligning sections of the query to the REF and ALT
        sequence extended a few bases to the left and right.

        If reference is None, alleles are detected by inspecting the
        existing alignment (via the CIGAR).

        chromosome -- name of chromosome to work on
        variants -- list of vcf.VcfVariant objects
        sample -- name of sample to work on. If None, read group information is
            ignored and all reads in the file are used.
        reference -- reference sequence of the given chromosome (or None)
        regions -- list of start,end tuples (end can be None)
        """
        # Since variants are identified by position, positions must be unique.
        if __debug__ and variants:
            varposc = Counter(variant.position for variant in variants)
            pos, count = varposc.most_common()[0]
            assert count == 1, "Position {} occurs more than once in variant list.".format(pos)

        logger.debug("Extracting Usable Alignments")
        alignments = self._usable_alignments(chromosome, sample)
        logger.debug("Converting Alignments to Read Objects")
        reads = self._alignments_to_reads(alignments, variants, sample, reference)
        logger.debug("Grouping Reads into ReadSet Object")
        grouped_reads = self._group_paired_reads(reads)
        logger.debug("ReadSet Object Successfully Created")
        readset = self._make_readset_from_grouped_reads(grouped_reads)
        return readset

    @staticmethod
    def _group_paired_reads(reads: Iterable[Read]) -> Iterator[List[Read]]:
        """
        Group reads into paired-end read pairs. Uses name, source_id and sample_id
        as grouping key.

        TODO
        Grouping by name should be sufficient since the SAM spec states:
        "Reads/segments having identical QNAME are regarded to come from the same template."
        """
        groups = defaultdict(list)
        for read in reads:
            groups[(read.source_id, read.name, read.sample_id)].append(read)
        for group in groups.values():
            if len(group) > 2:
                raise ReadSetError(
                    "Read name {!r} occurs more than twice in the input file".format(group[0].name)
                )
            yield group

    def _usable_alignments(self, chromosome, sample, regions=None):
        """
        Retrieve usable (suficient mapping quality, not secondary etc.)
        alignments from the alignment file
        """
        if regions is None:
            regions = [(0, None)]
        for s, e in regions:
            for alignment in self._reader.fetch(
                reference=chromosome, sample=sample, start=s, end=e
            ):
                # TODO handle additional alignments correctly!
                # find out why they are sometimes overlapping/redundant
                if (
                    alignment.bam_alignment.flag & 2048 != 0
                    or alignment.bam_alignment.mapping_quality < self._mapq_threshold
                    or alignment.bam_alignment.is_secondary
                    or alignment.bam_alignment.is_unmapped
                    or alignment.bam_alignment.is_duplicate
                ):
                    continue
                yield alignment

    def _alignments_to_reads(self, alignments, variants, sample, reference):
        """
        Convert BAM alignments to Read objects.

        If reference is not None, alleles are detected through re-alignment.

        Yield Read objects.
        """
        # FIXME hard-coded zero
        numeric_sample_id = 0 if sample is None else self._numeric_sample_ids[sample]
        if reference is not None:
            # Copy the pyfaidx.FastaRecord into a str for faster access
            reference = reference[:]
            normalized_variants = variants
        else:
            normalized_variants = variants

        i = 0  # index into variants
        for alignment in alignments:
            # Skip variants that are to the left of this read
            while (
                i < len(normalized_variants)
                and normalized_variants[i].position < alignment.bam_alignment.reference_start
            ):
                i += 1

            barcode = ""
            if alignment.bam_alignment.has_tag("BX"):
                barcode = alignment.bam_alignment.get_tag("BX")

            read = Read(
                alignment.bam_alignment.qname,
                alignment.bam_alignment.mapq,
                alignment.source_id,
                numeric_sample_id,
                alignment.bam_alignment.reference_start,
                barcode,
            )

            detected = self.detect_alleles_by_alignment(
                self._aligner,
                variants,
                i,
                alignment.bam_alignment,
                reference,
                self._realign_mode,
                self._overhang,
                self._em_params
            )
            for j, allele, em, quality in detected:
                read.add_variant(variants[j].position, allele, em, quality)
            if read:  # At least one variant covered and detected
                yield read

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        self._aligner.__dealloc__()
        self._reader.close()


def merge_two_reads(read1: Read, read2: Read) -> Read:
    """
    Merge two reads *that belong to the same haplotype* (such as the two
    ends of a paired-end read) into a single Read. Overlaps are allowed.
    """
    assert read1.is_sorted()
    assert read2.is_sorted()
    if read2:
        result = Read(
            read1.name,
            read1.mapqs[0],
            read1.source_id,
            read1.sample_id,
            read1.reference_start,
            read1.BX_tag,
        )
        result.add_mapq(read2.mapqs[0])
    else:
        return read1

    i1 = 0
    i2 = 0

    def add1():
        result.add_variant(read1[i1].position, read1[i1].allele, read1[i1].quality)

    def add2():
        result.add_variant(read2[i2].position, read2[i2].allele, read2[i2].quality)

    while i1 < len(read1) or i2 < len(read2):
        if i1 == len(read1):
            add2()
            i2 += 1
            continue
        if i2 == len(read2):
            add1()
            i1 += 1
            continue
        variant1 = read1[i1]
        variant2 = read2[i2]
        if variant2.position < variant1.position:
            add2()
            i2 += 1
        elif variant2.position > variant1.position:
            add1()
            i1 += 1
        else:
            # Variant on self-overlapping read pair
            assert read1[i1].position == read2[i2].position
            # If both alleles agree, merge into single variant and add up qualities
            if read1[i1].allele == read2[i2].allele:
                quality = read1[i1].quality + read2[i2].quality
                result.add_variant(read1[i1].position, read1[i1].allele, quality)
            else:
                # Otherwise, take variant with highest base quality and discard the other.
                if read1[i1].quality >= read2[i2].quality:
                    add1()
                else:
                    add2()
            i1 += 1
            i2 += 1
    return result


def merge_reads(*reads: Read) -> Read:
    """
    Merge multiple reads that belong to the same haplotype into a single Read.

    If the iterable is empty, a ValueError is raised.

    This 'naive' version just calls merge_two_reads repeatedly on all the reads.

    # TODO
    # The actual challenge is dealing with conflicts in variants covered by
    # more than one read. A solution would be to not merge if there are any
    # (or too many) conflicts and let the main algorithm deal with it.
    """
    it = iter(reads)
    try:
        read = next(it)
    except StopIteration:
        raise ValueError("no reads to merge")
    assert read.is_sorted()
    for partner in it:
        read = merge_two_reads(read, partner)
    return read
