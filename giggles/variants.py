"""
Detect variants in reads.
"""
# Code modified from WhatsHap (https://github.com/whatshap/whatshap)

import math
import re
from collections import defaultdict, Counter, namedtuple
from typing import Iterable, Iterator, List
from pywfa import WavefrontAligner

from giggles.logger import logger
from giggles.core import Read, ReadSet
from giggles.gaf import GafParser, rGFA, GafAlignment
from giggles.align import edit_distance
from giggles._variants import _iterate_cigar


class AlignmentReader:
    """
    Superclass for the GAF Readers (and any other alignment file readers)
    """
    def __init__(
            self,
            path: List[str],
            mapq_threshold: int,
            realign_mode: str,
            overhang: int,
            gap_start: int,
            gap_extend: int,
            default_mismatch: int):

        self._path = path
        self._mapq_threshold = mapq_threshold
        self._realign_mode = realign_mode
        if realign_mode == "edit":
            self._aligner = edit_distance
        elif realign_mode == "wfa":
            self._aligner = WavefrontAligner(mismatch=default_mismatch, 
                                         gap_opening=gap_start,
                                         gap_extension=gap_extend,
                                         scope='score')
        self._gap_start = gap_start
        self._gap_extend = gap_extend
        self._default_mismatch = default_mismatch
        self._overhang = overhang
        
    @property
    def n_paths(self):
        return len(self._paths)

    @staticmethod
    def _make_readset_from_grouped_reads(groups: Iterable[List[Read]], reg_const: int, base_const: float) -> ReadSet:
        read_set = ReadSet()
        for group in groups:
            if group is None:
                return None
            read_set.add(merge_reads(*group, reg_const = reg_const, base_const = base_const))
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
            read,
            cigartuples,
            i,
            consumed,
            query_pos,
            reference,
            mode,
            overhang):
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

        Return a tuple (allele, scores) where
        allele -- allele index with the max score
        scores -- list of alignment scores for all alleles (first is reference, then alternatives)
        """
        # Do not process symbolic alleles like <DEL>, <DUP>, etc.
        if any([alt.startswith("<") for alt in variant.alternative_allele]):
            return None, None

        # There is a big difference between the previous implementation and what is needed.
        # In the previous code, the CIGAR is against the reference always and hence we need to realign only for the alternate alleles.
        # With GAF, the CIGAR is not always against the reference (sometimes it is not against ref or any of the alt and can be with a path that is not an allele traversal)
        # So we need to generalize the process to realign using the variant record and the cigar tuples.
        left_cigar, right_cigar = AlignmentReader.split_cigar(cigartuples, i, consumed)

        if not variant.is_sv():
            # this is an external variant
            # overhang is set to 10
            left_ref_bases, left_query_bases = AlignmentReader.cigar_prefix_length(cigar=left_cigar[::-1], reference_bases=10)
            if variant.reference_allele == "*":
                ref_allele = ""
            else:
                ref_allele = variant.reference_allele
            # This should not be len(ref_allele)! This should be whatever path is followed in the node path!
            right_ref_bases, right_query_bases = AlignmentReader.cigar_prefix_length(cigar=right_cigar, reference_bases=variant.length_on_path + 10)

            assert variant.position - left_ref_bases >= 0
            assert variant.position + right_ref_bases <= len(reference)

            query = read.query_sequence[query_pos - left_query_bases : query_pos + right_query_bases]
            
            left_overhang = reference[variant.position - left_ref_bases : variant.position]
            right_overhang = reference[variant.position + right_ref_bases - 10 : variant.position + right_ref_bases]

            ref = left_overhang + ref_allele + right_overhang
            
            alts = []
            for alt_allele in variant.alternative_allele:
                if alt_allele != "*":
                    alt = left_overhang + alt_allele + right_overhang
                else:
                    alt = left_overhang + right_overhang
                alts.append(alt)

            scores = []
            max_score = -1e15
            max_allele = None
            for index, allele in enumerate([ref]+alts):
                scores.append(edit_distance(query, allele))
                if scores[index] > max_score:
                    max_score = scores[index]
                    max_allele = index

            return max_allele, scores

        # This is a SV variant
        left_ref_bases, left_query_bases = AlignmentReader.cigar_prefix_length(cigar=left_cigar[::-1], reference_bases=overhang)
        
        if variant.reference_allele == "*":
            ref_allele = ""
        else:
            ref_allele = variant.reference_allele

        # This should not be len(ref_allele)! This should be whatever path is followed in the node path!
        right_ref_bases, right_query_bases = AlignmentReader.cigar_prefix_length(cigar=right_cigar, reference_bases=variant.length_on_path + overhang)

        assert variant.position - left_ref_bases >= 0
        assert variant.position + right_ref_bases <= len(reference)

        query = read.query_sequence[query_pos - left_query_bases : query_pos + right_query_bases]
        
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

        scores = []
        max_score = -1e15
        max_allele = None
        for index, allele in enumerate([ref]+alts):
            if (abs(len(query) - len(allele)) > 5000 ) and (len(query)/len(allele) > 1.5 or len(query)/len(allele) < 1/1.5):
                # If the distance between the allele and query is too much, add a known low value
                scores.append(-1e10)
            else:
                if mode == "edit":
                    scores.append(aligner(query, allele))    #edit distance is positive.
                elif mode == "wfa":
                    scores.append(aligner(query, allele).score)
            if scores[index] > max_score:
                max_score = scores[index]
                max_allele = index

        return max_allele, scores

    @staticmethod
    def detect_alleles_by_alignment(
        aligner,
        variants,
        j,
        read,
        reference,
        mode,
        overhang=10,
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
            allele, scores = AlignmentReader.realign(
                aligner,
                variants[index],
                read,
                cigartuples,
                i,
                consumed,
                query_pos,
                reference,
                mode,
                overhang
            )

            if allele is not None:
                yield (index, allele, scores)


class GAFReader(AlignmentReader):
    """
    Associate VCF variant with GAF Read.
    """

    def __init__(
        self,
        alignment_files: List[str],
        reference: rGFA,
        read_fasta_files: List[str],
        mapq_threshold: int = 20,
        realign_mode: str = "wfa_full",
        overhang: int = 10,
        gap_start: int = 3,
        gap_extend: int = 1,
        default_mismatch: int = 2,
    ):
        super().__init__(
            alignment_files, 
            mapq_threshold, 
            realign_mode, 
            overhang, 
            gap_start, 
            gap_extend, 
            default_mismatch)

        self._reader = GafParser(alignment_files=alignment_files, reference=reference, read_fasta_files=read_fasta_files, mapq=self._mapq_threshold)

    def has_reference(self, chromosome):
        return self._reader.has_reference(chromosome)

    def read(self, chromosome, variants) -> ReadSet:
        """
        Detect alleles and return a ReadSet object containing reads representing
        the given variants.

        Using the provided reference, re-alignment is done which generates scores
        to be used in the HMM.

        chromosome -- name of chromosome to work on
        variants -- list of vcf.VcfVariant objects
        reference -- Here the variable does nothing. Kept to maintain compatibility with ReadSetReader
        """
        # Since variants are identified by position, positions must be unique.
        if __debug__ and variants:
            varposc = Counter(variant.position_on_ref for variant in variants)
            pos, count = varposc.most_common()[0]
            assert count == 1, "Position {} occurs more than once in variant list.".format(pos)

        logger.debug("Extracting Usable Alignments")
        alignments = self._usable_alignments(chromosome)
        logger.debug("Finding Variants in Alignments")
        updated_variants = self._update_variants_in_alignments(alignments, variants)
        logger.debug("Converting Alignments to Read Objects")
        reads = self._alignments_to_reads(updated_variants)
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
            if read is None:
                yield None
            if groups[(read.source_id, read.name)] == []:
                groups[(read.source_id, read.name)] = [read]       # Keeping this as a list so that I dont need to change _make_readset_from_grouped_reads()
            else:
                old_read = groups[(read.source_id, read.name)][0]

                # Check the number of variants it covers
                if len(old_read) < len(read):
                    groups[(read.source_id, read.name)] = [read]
                
                # Check the mapping quality
                if old_read.mapqs < read.mapqs:
                    groups[(read.source_id, read.name)] = [read]
        
        for group in groups.values():
            if len(group) > 1:
                raise Exception(f"Read name {group[0].name} occurs more than twice in the input file")
            yield group

    def _usable_alignments(self, chromosome):
        """
        Retrieve usable (suficient mapping quality, not secondary etc.)
        alignments from the alignment file
        """
        
        for alignment in self._reader(contig=chromosome):
            yield alignment

    @staticmethod
    def find_variants_in_alignment(alignment, variants, rgfa, variant_pointer):
        # finding the variants covered in the GAF alignment.
        
        '''
        These two functions use the path in the GAF alignment to find the start and end of the alignment on the reference.
        They also return the index of the first and last reference nodes in the path.

        What do I need to index for?
        - The main reason is to find whether the alignment covers partial bubbles in the start and end.
        - If the alignment starts or ends with a non-reference node, then it covers a partial bubble in the start or end.
        '''
        # function to get the start of the alignment on the reference
        # start_node_index returns the the index in alignment.path which is the first reference node in the path 
        alignment_start_on_ref, start_node_index, start_scaffold_node = GafAlignment.get_alignment_start_on_ref(alignment, rgfa)
        # function to get the end of the alignment on the reference
        # end_node_index returns the the index in alignment.path which is the first last node in the path 
        alignment_end_on_ref, end_node_index, end_scaffold_node = GafAlignment.get_alignment_end_on_ref(alignment, rgfa)

        logger.trace(f"Start of alignment on reference: {alignment_start_on_ref}")
        logger.trace(f"End of alignment on reference: {alignment_end_on_ref}")

        variants_in_alignment = []
        if alignment_start_on_ref is not None:
            # this checks if the end of the current variant is less than the start position of the next alignment
            while (variant_pointer+1 < len(variants) and (variants[variant_pointer].position_on_ref + len(variants[variant_pointer].reference_allele) < alignment_start_on_ref)):
                variant_pointer += 1
        else:
            node_in_path = alignment.path[1]  # first node in the path. does not matter which one since all should have same bo tag
            bo_tag = rgfa.get_node(node_in_path).tags['BO']
            if variant_pointer == 0:
                # need to check if the bubble-only alignments are before the first variant
                tmp_pointer = variant_pointer
                while tmp_pointer < len(variants) and variants[tmp_pointer].get_variant_bo(rgfa) is None:
                    tmp_pointer += 1
                if bo_tag == variants[tmp_pointer].get_variant_bo(rgfa):
                    # shifting variant pointer if the bubble-only alignment was for the first bubble variant
                    variant_pointer = tmp_pointer
            else:
                if variants[variant_pointer].get_variant_bo(rgfa) == None:
                    # variant pointer is currently at an external variant but we are looking at a alignment exclusively inside a bubble
                    variant_pointer += 1
                while variant_pointer + 1 < len(variants) and variants[variant_pointer].get_variant_bo(rgfa) < bo_tag:
                    variant_pointer += 1

        # the above definition of variant pointer cannot consider the case where the
        # first bubble is a partial bubble and the first node in the path is a non-
        # reference node.
        # Example 1: >sNR>sRS>sRNS where NR is non-ref, RS is ref-scaffold, and RNS is ref-non-scaffold.
        # in this case the start on alignment will become SO pos of sR-1 even though
        # the previous bubble is still covered.
        # the above issue does not happen in the case of >sNR>sRNS>sRS (Example 2)

        # this case will happen only when there is a scaffold node present
        if start_scaffold_node is not None:
            # if the start ref index is not the first node position
            if start_node_index > 1:
                # checking if the first ref node found is the scaffold node.
                # if yes, then this corresponds to Example 1. Variant pointer needs to be shifted back.
                # if no, then this corresponds to Example 2. No changes to variant pointer.
                if start_node_index == start_scaffold_node:
                    variant_pointer -= 1
                    if variant_pointer < 0:
                        # in case the variant pointer
                        variant_pointer = 0

        if alignment_start_on_ref is not None:
            assert alignment_end_on_ref is not None
            # This alignment ends before the first unprocessed variant starts
            if alignment_end_on_ref < variants[variant_pointer].position_on_ref:
                logger.trace(f'Alignment end: {alignment_end_on_ref}, < variants[variant_pointer].position_on_ref: {variants[variant_pointer].position_on_ref}')
                logger.trace(f'Alignment ends before the first unprocessed variant. No variants found in read {alignment.read_id}!')
                return None

            # This alignment starts after the last variant on the chromosome
            if (variant_pointer == len(variants) - 1) and (alignment_start_on_ref > variants[variant_pointer].position_on_ref + len(variants[variant_pointer].reference_allele)):
                logger.trace(f'Alignment starts after the last variant. No variants found in read {alignment.read_id}!\n')
                return None
        else:
            # If the alignment has no reference nodes and hence has no start or end on reference
            # now we need to look at the variants in terms of BO tags
            assert alignment_end_on_ref is None
            node_in_path = alignment.path[1]  # first node in the path. does not matter which one since all should have same bo tag
            bo_tag = rgfa.get_node(node_in_path).tags['BO']
            # the alignment ends before the first unprocessed variant
            if variants[variant_pointer].get_variant_bo(rgfa) is not None:
                # the variant pointer is currently at a SV variant
                if bo_tag < variants[variant_pointer].get_variant_bo(rgfa):
                    logger.trace(f'Alignment ends before the first unprocessed variant. No variants found in read {alignment.read_id}!')
                    return None
            else:
                # the variant pointer is currently at an external variant
                assert variant_pointer == 0, "This case should only happen when the alignment is before the first variant. There is some potential issues in sorting."
                logger.trace(f'Alignment ends before the first unprocessed variant. No variants found in read {alignment.read_id}!')
                return None
            # the alignment starts after the last variant in the chromosomes
            if (variant_pointer == len(variants) - 1) and (variants[variant_pointer].get_variant_bo(rgfa) == None):
                # if the last variant is a external variant. If this has been reached, then no variants can be found
                logger.trace(f'Alignment starts after the last variant. No variants found in read {alignment.read_id}!\n')
                return None
            if (variant_pointer == len(variants) - 1) and (bo_tag > variants[variant_pointer].get_variant_bo(rgfa)):
                # if the last variant is a bubble variant.
                logger.trace(f'Alignment starts after the last variant. No variants found in read {alignment.read_id}!\n')
                return None
            

        count_ext = 0   # Count of external variants
        count_sv = 0    # Count of structural variants
        end_pointer = variant_pointer
        variants_in_alignment.append(variants[variant_pointer])

        if start_scaffold_node is not None:
            # looking at the case where the alignment touches a scaffold node
            assert end_scaffold_node is not None
            # checking for the special case where the first variant has been added and the alignment has nodes before that.
            # this requires special attention since there is no bubble before the first variant
            partial_bubble_start_of_chrom = False
            if variant_pointer == 0 and start_scaffold_node > 1:
                # this condition alone cannot distinguish between the following two cases:
                # Case 1: Alignment actually starts at the beginning of the chromosome where the bubble is not defined.
                #         Let this case have an alignment like >sSB1>sSB2>sRS1>sRNS1>sRNS2>sRS2....
                #         As we can see, the variant pointer points to first variant and also start_scaffold_node (sRS1) is > 1 since it is the third node.
                # Case 2: Alignment starts in the middle of the first variant.
                #         Let this case have an alignment like >sRNS1>sRNS2>sRS2....
                #         As we can see, the variant pointer still points to first variant since it is being partially covered and the start_scaffold_node (sRS3) is > 1 since it is the third node.
                # Note: This case only happens when there are no external variants present in the first scaffold node.
                if not alignment_start_on_ref > variants[0].position_on_ref:
                    # to distinguish between the two cases, we look if the alignment starts before or after the first variant
                    partial_bubble_start_of_chrom = True

            if not variants_in_alignment[0].is_sv():
                variants[0].length_on_path = len(variants[0].reference_allele)
                count_ext += 1
            else:
                count_sv += 1
            while (end_pointer+1 < len(variants) and (variants[end_pointer+1].position_on_ref <= alignment_end_on_ref)):
                end_pointer += 1
                variants_in_alignment.append(variants[end_pointer])
                if not variants[end_pointer].is_sv():
                    variants[end_pointer].length_on_path = len(variants[end_pointer].reference_allele)  # This is the length of the variant on the alignment path.
                    count_ext += 1
                else:
                    count_sv += 1
            
            # checking for the special case where the last variant has been added and the alignment has nodes after that.
            partial_bubble_end_of_chrom = False
            if end_pointer == len(variants) - 1 and end_scaffold_node < len(alignment.path) - 1:
                # taking into consideration similar arguments for the end of the chromosome as given for the start
                # to distinguish between the two cases, we look if the alignment ends on reference after the reference allele of the last variant
                if not alignment_end_on_ref <= variants[-1].position_on_ref + len(variants[-1].reference_allele):
                    partial_bubble_end_of_chrom = True
        else:
            # looking at the case where the alignment does not touch any scaffold node
            assert end_scaffold_node is None
            partial_bubble_end_of_chrom = False
            partial_bubble_start_of_chrom = False
            assert variants[variant_pointer].is_sv(), "The variant has to be a bubble variant"
            count_sv += 1

        return variants_in_alignment, alignment_start_on_ref, variant_pointer, start_scaffold_node, end_scaffold_node, count_ext, count_sv, partial_bubble_start_of_chrom, partial_bubble_end_of_chrom

    def _update_variants_in_alignments(self, alignments, variants):
        """
        Finding the variants covered in the alignment and updating their positions on the path.

        Yields
        - variants_in_alignment: VcfVariant objects found in the alignment
        - alignment: GafAlignment object
        - alignment_start_on_ref: starting position of the alignment on reference
        - reference: reference seqeunce of the path given in the Gaf Alignment
        - partial_bubble_start_of_chrom_copy: if the start of the chromosome is present in the alignment (only needed for testing purposes)
        - partial_bubble_end_of_chrom: if the end of the chromosome is present in the alignment (only needed for testing purposes)
        """
        rgfa = self._reader._reference
        
        variant_pointer = 0     # Points to the variant which has not been processed in all Reads
        for alignment in alignments:
            # if no alignment found for the chromosome, return None
            if alignment is None:
                yield None
            # determine what variants are present in the alignment
            # returns an alignment in the correct orientation.
            alignment, _ = GafAlignment.check_reverse(alignment, rgfa)
            logger.debug(f"Processing alignment {alignment.read_id} on {alignment.source_id}")
            logger.trace(f"Alignment Path: {''.join(alignment.path)}")
            result = GAFReader.find_variants_in_alignment(alignment, variants, rgfa, variant_pointer)
            # in the case None is returned.
            if result is None:
                continue
            # extract items from result
            variants_in_alignment, alignment_start_on_ref, variant_pointer, start_scaffold_node_index_on_path, end_scaffold_node_index_on_path, count_ext, count_sv, partial_bubble_start_of_chrom, partial_bubble_end_of_chrom = result
            # need a copy of this since we are setting partial_bubble_start_of_chrom to False if it is True.
            partial_bubble_start_of_chrom_copy = partial_bubble_start_of_chrom
            logger.trace(f"Found variants in alignment: {variants_in_alignment}")
            logger.trace(f"Number of variants: {len(variants_in_alignment)}, #EXT: {count_ext}, #SV: {count_sv}")
            logger.trace(f"Detected partial bubbles at chromosome ends: start={partial_bubble_start_of_chrom}, end={partial_bubble_end_of_chrom}")
            assert (count_sv + count_ext) == len(variants_in_alignment), "The number of SV and EXT bubbles in the alignment does not match the number of variants found in the alignment."

            # initializing booleans for path starts and ends with scaffold
            path_starts_with_scaffold = False
            path_ends_with_scaffold = False
            if start_scaffold_node_index_on_path == 1:
                path_starts_with_scaffold = True
                # This implies that there are no partial bubbles at the start of the alignment
            if end_scaffold_node_index_on_path == len(alignment.path) - 1:
                path_ends_with_scaffold = True
                # This implies that there are no partial bubbles at the end of the alignment

            # adding the lengths of the variants on the alignment path for SV bubble variants
            reference = ""
            len_on_path = 0     # This keeps track of the length of the allele in the alignment path for different bubbles
            active_pointer = 0  # This keeps track of the variant in the variants_in_alignment list.
            sv_counter = 0      # counting sv bubbles (for comparing later to the results of find_variants_in_alignment)
            ext_counter = 0     # counting ext bubbles (for comparing later to the results of find_variants_in_alignment)
            for index, n in enumerate(alignment.path):
                if n in ['>', '<']:
                    orient = n
                    continue
                node = rgfa.get_node(n)
                # keeping track of reference length before updating (to identify where the bubbles start)
                prev_ref_len = len(reference)
                # Reference updated
                if orient == '<':
                    reference += GAFReader.reverse_complement(node.sequence)
                else:
                    reference += node.sequence
                while (active_pointer < len(variants_in_alignment)) and (not variants_in_alignment[active_pointer].is_sv()):
                    # if the current variant is an extension variant
                    # going to the next variant since we only need the SV bubble variants
                    ext_counter += 1
                    # logger.debug(f"Found EXT variant {variants_in_alignment[active_pointer].id}. Now active_pointer is {active_pointer+1}")
                    active_pointer += 1
                if index == 1 and path_starts_with_scaffold:
                    # this is the case where the path starts with a scaffold
                    # hence the first node needs to be skipped otherwise the first bubble will be tagged with length_on_path = 0.
                    assert node.tags['NO'] == 0, "Path starts with scaffold but the first node is not a scaffold node."
                    continue
                # this is a non-scaffold node
                if node.tags['NO'] != 0:
                    len_on_path += len(node.sequence)
                # this is a scaffold node
                if node.tags['NO'] == 0:
                    if partial_bubble_start_of_chrom:
                        # the alignment spans the start of the chromosome where the bubble is not defined.
                        logger.trace(f"Partial bubble found at the start of the chromosome. Bubble variant is not defined. Skipping.")
                        len_on_path = 0
                        partial_bubble_start_of_chrom = False   # setting it as false
                        continue
                    if index == start_scaffold_node_index_on_path and index != 1:
                        # this is the case where the alignment spans a partial bubble in the beginning
                        # in this case, the start of the alignment on the path has to be taken into consideration
                        assert active_pointer == 0, "Active pointer should be at the start when considering partial bubble."
                        variants_in_alignment[active_pointer].length_on_path = len_on_path - alignment.p_start
                        variants_in_alignment[active_pointer].position = alignment.p_start
                        active_pointer += 1
                        sv_counter += 1
                        len_on_path = 0
                        continue
                    # this is a scaffold node that is not at the end of the alignment
                    # this implies a bubble has ended and the next one will begin.
                    # the +1 in the length is to keep in consideration that the bubble variants always contain the last nucleotide of the previous node. Done to avoid blank alleles when there are no nodes between the scaffold node.
                    # the -1 in the position is to consider the same and shift the start of the allele by 1bp
                    variants_in_alignment[active_pointer].length_on_path = len_on_path + 1
                    variants_in_alignment[active_pointer].position = prev_ref_len - len_on_path - 1     # the position where the variant starts is the point where the reference (before being updated in this node step) ends minus the length the variant has on the path
                    len_on_path = 0
                    #logger.debug(f"Found SV variant {variants_in_alignment[active_pointer].id}. Now active_pointer is {active_pointer+1}")
                    active_pointer += 1
                    sv_counter += 1
                    continue
                if (not path_ends_with_scaffold) and (index == len(alignment.path) - 1):
                    if partial_bubble_end_of_chrom:
                        # this is the case where there is a partial alignment to the final bubble on the chromosome
                        # but the bubble is not defined due to no scaffold nodes at the end.
                        logger.trace(f"Partial bubble found at the end of the chromosome. Bubble variant is not defined. Skipping.")
                        assert active_pointer == len(variants_in_alignment), "Active pointer should be out of bounds since there are no more variants"
                        continue
                    # if the alignment does not end with a scaffold node and has a partial alignment to a defined bubble.
                    # the length on that path still needs to be stored.
                    # this also requires to take into consideration where the alignment ends on the path
                    assert active_pointer == len(variants_in_alignment)-1, "Active pointer is out of bounds when storing alignment to partial bubble at the end."
                    assert sv_counter == count_sv - 1, "Active pointer is not at the last SV when storing alignment to partial bubble at the end."
                    logger.trace(f"Found SV variant {variants_in_alignment[active_pointer].id} as partial bubble alignment.")
                    if start_scaffold_node_index_on_path is not None:
                        # the +1 in the length is to keep in consideration that the bubble variants always contain the last nucleotide of the previous node. Done to avoid blank alleles when there are no nodes between the scaffold node.
                        # the -1 in the position is to consider the same and shift the start of the allele by 1bp
                        variants_in_alignment[active_pointer].length_on_path = len_on_path - (alignment.p_len - alignment.p_end) + 1
                        variants_in_alignment[active_pointer].position = len(reference) - len_on_path - 1
                    else:
                        # this is the case where the alignment is only inside a bubble
                        assert len(variants_in_alignment) == 1, "There should be only one variant in the alignment since it is only inside a bubble"
                        variants_in_alignment[active_pointer].length_on_path = alignment.p_end - alignment.p_start + 1
                        variants_in_alignment[active_pointer].position = alignment.p_start
                    sv_counter += 1
            # only asserting for sv counter.
            # does not make sense for ext variants since the ext variants after the last sv are not considered.
            assert sv_counter == count_sv, "The number of SV bubbles in the alignment does not match the number of SV bubbles found in the alignment."

            # determining the position of each variant in the path
            # if the alignment has partial bubble alignment (not at the start of the chromosome)
            position_tracker = None
            for index, variant in enumerate(variants_in_alignment):
                # this is an SV variant. its information should already be stored.
                if variant.is_sv():
                    position_tracker = variant.position + variant.length_on_path
                    continue
                if index == 0:
                    # the first variant
                    assert not variant.is_sv(), "The first variant in the alignment cannot be a SV variant"
                    # determine external variants position in its scaffold node
                    start_scaffold_node = alignment.path[start_scaffold_node_index_on_path]
                    assert ((variant.position_on_ref >= rgfa.get_node(start_scaffold_node).start) and (variant.position_on_ref < rgfa.get_node(start_scaffold_node).start+rgfa.get_node(start_scaffold_node).tags['LN'])), "First external variant in the alignment is found out of bounds from the first scaffold node"
                    position_on_scaffold = variant.position_on_ref - rgfa.get_node(start_scaffold_node).start
                    position_tracker = 0
                    # iterating through all the nodes left of the start scaffold node and getting their path sum.
                    for node in alignment.path[:start_scaffold_node_index_on_path]:
                        if node in ['>', '<']:
                            continue
                        position_tracker += rgfa.get_node(node).tags['LN']
                    position_tracker += position_on_scaffold
                    variant.position = position_tracker
                    continue
                if not variants_in_alignment[index-1].is_sv():
                    # if the previous variant was an external variant
                    # the position of the current external variant simply is the distance between the previous the last ext variant and this one
                    position_tracker += variant.position_on_ref - variants_in_alignment[index-1].position_on_ref
                    variant.position = position_tracker
                else:
                    # the external variant should be on the scaffold node given in the id of the previous SV variant
                    scaffold_node = variants_in_alignment[index-1].id.split('>')[-1]
                    assert ((variant.position_on_ref >= rgfa.get_node(scaffold_node).start) and (variant.position_on_ref < rgfa.get_node(scaffold_node).start+rgfa.get_node(scaffold_node).tags['LN'])), "External variant in the alignment is found out of bounds from the scaffold node given in previous SV id"
                    # updating the position tracker with the position of the variant on the scaffold node
                    position_tracker += variant.position_on_ref - rgfa.get_node(scaffold_node).start
                    variant.position = position_tracker
            logger.trace(f"Positions of Variants in Alignment: {', '.join(f'{variant.id}: {variant.position}' for variant in variants_in_alignment)}")

            #if alignment.read_id == 'read2':
            #    exit()

            yield (variants_in_alignment, alignment, alignment_start_on_ref, reference, partial_bubble_start_of_chrom_copy, partial_bubble_end_of_chrom)

    def _alignments_to_reads(self, updated_variants):
        """
        Convert GAF alignments to Read objects.

        Yields Read objects
        """
        cg_letter_to_op = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, 'X': 7, '=': 8}
        Alignment = namedtuple('Alignment', ['cigartuples', 'reference_start', 'query_sequence'])        # Class created to maintain compatibility with old code
        
        for result in updated_variants:

            # if no alignment was found 
            if result is None:
                yield None

            variants_in_alignment, alignment, alignment_start_on_ref, reference, _, _ = result

            # Extract the aligned segement from the complete read sequence and create a new object.
            # Need cigartuples, and reference_start (where it starts in the reference. So the path start in this case.)
            #print(f'num_var_alignment: {len(variants_in_alignment)}')
            gaf_aligned_segment = alignment.sequence[alignment.q_start:alignment.q_end]
            cg_tuples = []
            cg = list(filter(None, re.split("([MIDNSHP=X])", alignment.cigar)))
            for i in range(0,len(cg),2):
                l = int(cg[i])
                op = cg_letter_to_op[cg[i+1]]
                cg_tuples.append((op,l))
            
            # This new variable is created to make the gaf alignments compatible with the old code.
            processed_alignment = Alignment(cigartuples=cg_tuples, reference_start=alignment.p_start, query_sequence=gaf_aligned_segment)
            
            read = Read(
                alignment.read_id,
                alignment.mapping_quality,
                alignment.source_id,
                alignment_start_on_ref,
            )
            
            detected = self.detect_alleles_by_alignment(
                self._aligner,
                variants_in_alignment,
                0,                              # Here this has been hardcoded to 0. In original code, this was the index of the first variant (index in the big list of variants) in the read. But now we have new list of variants just for this alignment.
                processed_alignment,
                reference,
                self._realign_mode,
                self._overhang)
            for j, allele, scores in detected:
                read.add_variant(variants_in_alignment[j].position_on_ref, allele, scores)
            if read:  # At least one variant covered and detected
                yield read

    def __enter__(self):
        return self

    def __exit__(self, *args):
        logger.debug("Closing GAFReader")
        self.close()

    def close(self):
        if type(self._aligner) is WavefrontAligner:
            logger.debug("Deallocating WavefrontAligner")
            self._aligner.__dealloc__()
        logger.debug("Closing GAFParser")
        self._reader.close()

# TODO: Do I need this?
def merge_two_reads(read1: Read, read2: Read, reg_const, base_const) -> Read:
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
            read1.reference_start,
            read1.BX_tag,
            reg_const,
            base_const
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

# TODO: Do I need this?
def merge_reads(*reads: Read, reg_const, base_const) -> Read:
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
        read = merge_two_reads(read, partner, reg_const, base_const)
    return read
