"""
Detect variants in reads.
"""
# Code modified from WhatsHap (https://github.com/whatshap/whatshap)

import re
from collections import defaultdict, Counter, namedtuple
from typing import Iterable, Iterator, List

from giggles.logger import logger
from giggles.core import Read, ReadSet
from giggles.gaf import GafParser, rGFA, GafAlignment
from giggles.vcf import VcfVariant
from giggles.align import edit_distance
from giggles.ext import WFAWrapper
from giggles._variants import _iterate_cigar


class Realigner:
    """
    Class to encapsulate different aligners
    """
    def __init__(self, bandwidth: int):
        self.bandwidth = bandwidth
        self.realigner = WFAWrapper(bandwidth)


    def get_distance(self, query: str, allele: str, state: int):
        """Calculates the distance between query and padded allele.

        Args
            query: the sequence on the read
            allele: the sequence of alleles given in the VCF (padded with some overhangs)
            state: type of alignment (0, 1, 2, or 3)
        
        Returns
            int: the distance between the two string (always positive)
        """
        dist = self.realigner.align(query, allele, state=state)
        return min(dist, self.bandwidth)
            

    def close(self):
        self.realigner.__dealloc__()



class AlignmentReader:
    """
    Superclass for the GAF Readers (and any other alignment file readers)
    """
    def __init__(
            self,
            path: List[str],
            mapq_threshold: int,
            bandwidth:int,
            overhang: int):

        self._path = path
        self._mapq_threshold = mapq_threshold
        self._aligner = Realigner(bandwidth)
        self._bandwidth = bandwidth
        self._overhang = overhang
        
    @property
    def n_paths(self):
        return len(self._paths)

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
        is_custom_graph: bool = False,
        bandwidth: int = 30,
        overhang: int = 10
    ):
        super().__init__(
            alignment_files, 
            mapq_threshold, 
            bandwidth,
            overhang)

        self._is_custom_graph = is_custom_graph
        self._reader = GafParser(alignment_files=alignment_files, reference=reference, read_fasta_files=read_fasta_files, mapq=self._mapq_threshold)

    def has_reference(self, chromosome):
        return self._reader.has_reference(chromosome)

    def read(self, chromosome: str, variants: List[VcfVariant]) -> ReadSet:
        """Detect alleles and return a ReadSet object containing reads representing the given variants.

        Args:
            chromosome: name of the chromosome
            variants: the variants present on the chromosome

        Returns:
            ReadSet: the ReadSet object containing the information about the variants found in each
                alignment and their distance scores from available alleles.
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
        readset = ReadSet()
        for group in grouped_reads:
            if group is None:
                continue
            readset.add(group[0])
        logger.debug("ReadSet Object Successfully Created")
        return readset      

    @staticmethod
    def _remove_duplicate_reads(reads: Iterable[Read]) -> Iterator[List[Read]]:
        """Removes reads which have been mapped multiple times and selects one best read.

        Args:
            reads: A list of Read objects

        Yields:
            dict{(int, str): Read}: a dictionary containing:
                - key is a tuple of the source id (which is an integer index of the file in which the read is present) and name of the read.
                - the value is the Read object created in _alignments_to_reads
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

    def _usable_alignments(self, chromosome: str) -> Iterator[GafAlignment]:
        """"Retrieves usable alignments from the alignment file.

        Currently does not have restrictions. In the future, this is where alignment 
        filtering will happen.

        Args:
            chromosome: The name of the chromosome.

        Yields:
            GafAlignment: The next usable alignment object.
        """
        
        for alignment in self._reader(contig=chromosome):
            yield alignment

    @staticmethod
    def find_variants_in_alignment(alignment: GafAlignment, variants: List[VcfVariant], rgfa: rGFA, variant_pointer: int):
        """Identifies variants covered by a specific alignment and assigns coverage states.

        This method synchronizes the alignment's genomic interval with a sorted list of 
        variants. It determines if the alignment fully covers, partially covers (enters/exits), 
        or is contained within variant bubbles (given by the VcfVariant.state variable).

        Args:
            alignment: The alignment object to process.
            variants: A sorted list of variant objects (0-based coordinates).
            rgfa: The reference graph object.
            variant_pointer: The index in the `variants` list to start searching from 
                (optimization for sorted data).

        Returns:
            tuple: A tuple containing 7 elements:
                1. list[VcfVariant]: The list of variants found in this alignment.
                2. int: Alignment start on reference (0-based, inclusive).
                3. int: Updated variant_pointer for the next iteration.
                4. int: Index of the first scaffold node in the alignment path.
                5. int: Index of the last scaffold node in the alignment path.
                6. int: Count of External (linear) variants found.
                7. int: Count of Structural (bubble) variants found.

            Returns a tuple of (None, ...) values if the pointer is out of bounds or 
            no overlap is possible.
        """
        
        # ------------------------------------------------------------------
        # STEP 1: Determine Alignment Coordinates & Topology
        # ------------------------------------------------------------------
        # Coordinates are 0-based, Half-Open [Start, End)
        alignment_start_on_ref, start_node_index, start_scaffold_node = GafAlignment.get_alignment_start_on_ref(alignment, rgfa)
        alignment_end_on_ref, end_node_index, end_scaffold_node = GafAlignment.get_alignment_end_on_ref(alignment, rgfa)

        '''
        The code below detects the following structure
        >sNR>sRS>.......>sRS>sNR
        --------        --------
            |               |
            v               v
        Non-reference nodes followed/preceded by reference scaffold nodes at the ends of the alignment
        sNR -> non-reference nodes (Note that this is specific to non-reference nodes and not reference non-scaffold nodes)
        sRS -> reference scaffold nodes
        '''
        has_non_ref_start = False
        if start_node_index is not None:
            if (start_node_index > 1) and (start_node_index == start_scaffold_node):
                has_non_ref_start = True
        has_non_ref_end = False
        if end_node_index is not None:
            if (end_node_index < len(alignment.path) - 1) and (end_node_index == end_scaffold_node):
                has_non_ref_end = True
        
        '''
        The code below detects the following structure
        >sNS>sRS>.......>sRS>sNS
        --------        --------
            |               |
            v               v
        Non-scaffold nodes followed/preceded by reference scaffold nodes at the ends of the alignment
        sNS -> non-scaffold nodes (Note that this is non-specific. sNS can be non-reference nodes and also reference non-scaffold nodes)
        sRS -> reference scaffold nodes
        '''
        has_partial_start = False
        if start_scaffold_node is not None and start_scaffold_node > 1:
            has_partial_start = True
        has_partial_end = False
        if end_scaffold_node is not None and end_scaffold_node < len(alignment.path) - 1:
            has_partial_end = True

        logger.trace(f"Start on ref: {alignment_start_on_ref}, End on ref: {alignment_end_on_ref}")

        # ------------------------------------------------------------------
        # STEP 2: Fast-Forward Pointer
        # ------------------------------------------------------------------
        if alignment_start_on_ref is not None:
            while variant_pointer < len(variants):
                variant = variants[variant_pointer]
                # Calculate Variant End (Exclusive)
                # Assumes variant.position_on_ref is already 0-based/anchor-free
                var_end_exclusive = variant.position_on_ref + len(variant.reference_allele)
                
                # Case 1: Variant is completely behind the alignment start
                if var_end_exclusive < alignment_start_on_ref:
                    variant_pointer += 1
                    continue
                
                # Case 2: Boundary Touch (Variant End == Alignment Start)
                if var_end_exclusive == alignment_start_on_ref:
                    if has_non_ref_start:
                        # We came from the bubble associated with this variant. Keep it.
                        break 
                    else:
                        # We started cleanly on the reference after this variant. Skip it.
                        variant_pointer += 1
                        continue
                
                # Case 3: Overlap (var_end > alignment_start)
                break
        else:
            # Bubble-Only Logic: Match via BO (Bubble Origin) tags
            node_in_path = alignment.path[1]
            alignment_bo = rgfa.get_node(node_in_path).tags.get('BO')
            
            while variant_pointer < len(variants):
                variant_bo = variants[variant_pointer].get_variant_bo(rgfa)
                if variant_bo is None or (alignment_bo is not None and variant_bo < alignment_bo):
                    # variant_bo is None: external variants
                    # variant_bo < alignment_bo: SV variants before the detected alignment
                    variant_pointer += 1
                else:
                    assert variant_bo == alignment_bo
                    break



        # ------------------------------------------------------------------
        # STEP 3: Safety Check (Overlap Validation)
        # ------------------------------------------------------------------
        empty_result = (None, alignment_start_on_ref, variant_pointer, start_scaffold_node, end_scaffold_node, 0, 0)
        
        if variant_pointer >= len(variants):
             return empty_result

        first_variant = variants[variant_pointer]
        
        if alignment_start_on_ref is not None:
            # Check if alignment ends before variant starts
            if alignment_end_on_ref < first_variant.position_on_ref:
                logger.trace(f'Alignment ends at {alignment_end_on_ref} before variant starts at {first_variant.position_on_ref}')
                return empty_result
            
            # Check Boundary: Alignment ends exactly where variant starts
            if alignment_end_on_ref == first_variant.position_on_ref and not has_non_ref_end:
                logger.trace(f'Alignment ends at {alignment_end_on_ref} before variant starts at {first_variant.position_on_ref} (Alignment end is index-exclusive)')
                return empty_result
        else:
            # Bubble-Only Validation
            assert alignment_end_on_ref is None
            alignment_bo = rgfa.get_node(alignment.path[1]).tags.get('BO')
            first_variant_bo = first_variant.get_variant_bo(rgfa)
            if first_variant_bo is None:
                # Pointer is at a linear variant, but alignment is in a bubble. No overlap.
                # (Removed the dangerous 'assert variant_pointer == 0' here)
                return empty_result
                
            if alignment_bo != first_variant_bo:
                logger.trace(f'BO Tag Mismatch: Alignment {alignment_bo} vs Variant {first_variant_bo}')
                return empty_result
        
        # ------------------------------------------------------------------
        # STEP 4: Collect Variants & Assign States
        # ------------------------------------------------------------------
        variants_in_alignment = []
        count_ext = 0
        count_sv = 0
        
        # Special Case: Bubble-Only Alignment (State 3)
        if start_scaffold_node is None:
            variants_in_alignment.append(first_variant)
            first_variant.state = 3
            assert first_variant.is_sv(), "Bubble-only alignment must map to an SV"
            count_sv += 1
            return variants_in_alignment, alignment_start_on_ref, variant_pointer, start_scaffold_node, end_scaffold_node, count_ext, count_sv

        # General Case: Iterate and Collect
        curr_idx = variant_pointer
        while curr_idx < len(variants):
            candidate = variants[curr_idx]
            cand_start = candidate.position_on_ref
            
            # Stop if candidate starts strictly after alignment ends
            if cand_start > alignment_end_on_ref:
                break
            
            # Stop if candidate starts EXACTLY at alignment end...
            if cand_start == alignment_end_on_ref:
                # ...UNLESS we exit via the bubble (has_non_ref_end)
                if has_non_ref_end and candidate.is_sv():
                    variants_in_alignment.append(candidate)
                    break # This is the last one
                else:
                    break # Boundary touched, but not entered
            
            # Standard Overlap
            variants_in_alignment.append(candidate)
            curr_idx += 1

        # Assign States
        for i, variant in enumerate(variants_in_alignment):
            if not variant.is_sv():
                variant.state = 0
                variant.length_on_path = len(variant.reference_allele)
                count_ext += 1
            else:
                count_sv += 1
                variant.state = 0 # Default: Full Coverage
                
                is_first = (i == 0)
                is_last = (i == len(variants_in_alignment) - 1)
                
                # State 1: Partial Start (Entered bubble from left)
                if is_first and has_partial_start:
                    variant.state = 1
                
                # State 2: Partial End (Exited bubble to right)
                if is_last and has_partial_end:
                    variant.state = 2
                    
        return variants_in_alignment, alignment_start_on_ref, variant_pointer, start_scaffold_node, end_scaffold_node, count_ext, count_sv


    def _calculate_sv_attributes(self, alignment: GafAlignment, variants_in_alignment: List[VcfVariant], rgfa: rGFA, start_scaf_idx: int):
        """Constructs the alignment sequence and calculates attributes for SV bubbles.

        Iterates through the alignment path nodes to build the full sequence string. 
        Simultaneously measures the physical length of Structural Variants (bubbles) 
        as they appear on the read path and updates their attributes in-place.

        Args:
            alignment: The alignment object.
            variants_in_alignment: The list of variants associated 
                with this alignment.
            rgfa: The reference graph object.
            start_scaf_idx: The index of the first scaffold node in the path 
                (used to handle partial start bubbles).

        Returns:
            str: The constructed DNA sequence of the alignment path.

        Note:
            This function modifies the `VcfVariant` objects within `variants_in_alignment` 
            in-place, setting their `length_on_path` and `position` attributes.
        """

        logger.trace(f'Setting attributes for SV bubble variants on {alignment.read_id}')
        reference_seq = ""
        len_on_path = 0     # Accumulator for current bubble length
        active_pointer = 0  # Pointer to variants_in_alignment list
        
        # Helper to skip EXT variants in the list (we only measure SVs here)
        def advance_pointer_to_next_sv():
            nonlocal active_pointer
            while active_pointer < len(variants_in_alignment) and not variants_in_alignment[active_pointer].is_sv():
                active_pointer += 1

        advance_pointer_to_next_sv()
        logger.trace(f'Finding the first SV bubble in the alignment. Active pointer is at {active_pointer}.')

        for index, node_id in enumerate(alignment.path):
            if node_id in ['>', '<']:
                orient = node_id
                continue

            logger.trace(f'Processing node {node_id} at index {index} of alignment.')
            
            node = rgfa.get_node(node_id)
            node_seq = node.sequence
            if orient == '<':
                node_seq = GAFReader.reverse_complement(node_seq)
            
            # Track start of this node in the built sequence
            current_seq_pos = len(reference_seq)
            reference_seq += node_seq
            
            is_scaffold = (node.tags.get('NO') == 0)

            # --- CASE A: Non-Scaffold Node (Inside Bubble) ---
            if not is_scaffold:
                logger.trace(f'Node {node_id} is a non-scaffold node. Adding its length to the accumulator of current bubble length.')
                len_on_path += len(node_seq)
                continue

            # --- CASE B: Scaffold Node (Bubble Boundary) ---
            
            # 1. Handle Start of Alignment (First Scaffold)
            # If we started with a partial bubble, we process it NOW.
            if index == start_scaf_idx:
                # Logic: We are at the first anchor. If there was a bubble before us,
                # it was a Partial Start (State 1).
                # Note: If start_scaf_idx is 1 (meaning path is >Ref...), len_on_path is 0. 
                # If start_scaf_idx > 1 (meaning path is >Bub>Ref...), len_on_path > 0.
                
                if len_on_path > 0 or (active_pointer < len(variants_in_alignment) and variants_in_alignment[active_pointer].state == 1):
                    if active_pointer < len(variants_in_alignment):
                        sv = variants_in_alignment[active_pointer]

                        logger.trace(f'Setting attributes of the first SV bubble ({sv}).')
                        # Partial Length: Total calculated - start_offset
                        sv.length_on_path = len_on_path - alignment.p_start
                        sv.position = alignment.p_start # Starts at beginning of read
                        
                        active_pointer += 1
                        advance_pointer_to_next_sv()
                        logger.trace(f'Finding the next SV bubble. Active pointer is at {active_pointer}.')
                        len_on_path = 0
                else:
                    # Clean start on scaffold, reset accumulator just in case
                    len_on_path = 0
                continue

            # 2. Handle Normal Bubble Closure
            # We hit a scaffold node, and it's not the start. 
            # This closes any bubble accumulating before this node.
            
            if active_pointer < len(variants_in_alignment):
                sv = variants_in_alignment[active_pointer]
                logger.trace(f'Setting attributes of the last SV bubble ({sv}) with complete coverage.')
                # Length: The sum of non-ref nodes we just traversed
                sv.length_on_path = len_on_path 
                
                # Position: Start of current node minus the bubble length
                # (This points to the index in reference_seq where the bubble began)
                sv.position = current_seq_pos - len_on_path
                
                active_pointer += 1
                advance_pointer_to_next_sv()
                logger.trace(f'Finding the next SV bubble. Active pointer is at {active_pointer}.')
            
            len_on_path = 0

        # --- CASE C: End of Alignment (Partial End) ---
        # If we finished the loop and still have len_on_path, or we are in State 2/3
        if active_pointer < len(variants_in_alignment):
            logger.trace(f'[3] End of Alignment. Active pointer is at {active_pointer}')
            sv = variants_in_alignment[active_pointer]
            logger.trace(f'Setting attributes of the last SV bubble ({sv}) which is covered partially.')
            
            # SAFETY CHECK: 
            # If this is State 0, it SHOULD have been closed by a scaffold node inside the loop.
            # If we are here, something is wrong with the State assignment or the Path logic.
            assert sv.state != 0, f"State 0 variant {sv.id} was not processed inside the loop! Path may be malformed."
            
            # Verify we are at the last SV
            # Logic: Length is whatever we accumulated, minus the unused end clip
            # formula: len_on_path - (path_len - p_end)
            
            distance_from_end = alignment.p_len - alignment.p_end
            
            # Special Handling for State 3 (Bubble Only)
            if sv.state == 3:
                 sv.length_on_path = alignment.p_end - alignment.p_start
                 sv.position = alignment.p_start
            else:
                # State 2 (Partial End)
                sv.length_on_path = len_on_path - distance_from_end
                sv.position = len(reference_seq) - len_on_path
            
        return reference_seq


    def _interpolate_ext_positions(self, variants: List[VcfVariant], alignment: GafAlignment, rgfa: rGFA, start_scaf_idx: int):
        """Calculates and updates positions for external (linear) variants on the read path.

        Since external variants (SNPs) do not correspond to graph nodes, their positions 
        are calculated by interpolating from the nearest anchor (either the start of the 
        alignment or the end of the previous Structural Variant).

        Args:
            variants: The list of variants in the alignment.
            alignment: The alignment object.
            rgfa: The reference graph object.
            start_scaf_idx: The index of the first scaffold node in the path.

        Returns:
            None: This function modifies the variant objects in the input list in-place.
        """

        logger.trace(f'Interpolating the positions of the external variants covered by {alignment.read_id}')
        current_path_pos = 0
        
        if start_scaf_idx is None:
            # if there are no scaffold nodes in the alignment,
            # then there are no external variants
            logger.trace('No scaffold nodes exist in this alignment. No external variants should be found.')
            return

        # Determine the initial anchor (Start Node)
        # Calculate path length up to the first scaffold node
        # (This handles the prefix if the read started in a bubble)
        prefix_len = 0
        for i in range(1, start_scaf_idx, 2):
            # Simple lookup, assuming valid path structure
            n_id = alignment.path[i]
            prefix_len += rgfa.get_node(n_id).tags['LN']
        current_path_pos = prefix_len
    
        # We need the start node object to calculate offsets
        start_node_id = alignment.path[start_scaf_idx]
        start_node = rgfa.get_node(start_node_id)
        start_node_start_ref = start_node.start

        for i, variant in enumerate(variants):
            if variant.is_sv():
                # Update anchor: The end of this SV becomes the new anchor
                current_path_pos = variant.position + variant.length_on_path
                continue
            
            # It's an External Variant
            if i == 0 or (i > 0 and variants[i-1].is_sv() and variants[i-1].state != 3):
                # CASE 1: Anchored to Start Node (First variant, or after an SV)
                # Wait, if i > 0, we should anchor to the previous SV's ID tag if possible?
                # The original code used `start_scaffold_node` for the first one.
                
                if i == 0:
                    # Offset on the scaffold node
                    offset = variant.position_on_ref - start_node_start_ref
                    assert current_path_pos == 0    # since this is the first variant and is external, then there should not be any cumulative path length stored
                    variant.position = offset
                else:
                    # Logic from original: Use previous SV ID to find scaffold
                    # ID format >s1>s2 -> Last one is scaffold
                    prev_sv = variants[i-1]
                    scaffold_id = prev_sv.id.split('>')[-1] 
                    scaffold = rgfa.get_node(scaffold_id)
                    
                    offset = variant.position_on_ref - scaffold.start
                    variant.position = current_path_pos + offset
            else:
                # CASE 2: Consecutive Ext Variants
                # Just add the delta from the previous Ext variant
                prev_var = variants[i-1]
                delta = variant.position_on_ref - prev_var.position_on_ref
                variant.position = prev_var.position + delta



    def _update_variants_in_alignments(self, alignments: Iterator[GafAlignment], variants: List[VcfVariant]):
        """Processes a stream of alignments to find, map, and attribute variants.

        Acts as the main orchestrator loop. It maintains a global pointer to the 
        sorted variant list to efficiently process sequential alignments.

        Args:
            alignments: An iterator yielding alignment objects.
            variants: A sorted list of variants for the current chromosome.

        Yields:
            tuple: A tuple containing:
                - list[VcfVariant]: The variants found in the alignment (with updated attributes).
                - GafAlignment: The processed alignment object.
                - int: The alignment's start position on the reference.
                - str: The reconstructed sequence of the alignment path.

            Yields `None` if the input alignment is None.
        """
        
        rgfa = self._reader._reference
        variant_pointer = 0  # Global pointer for the sorted variants list

        for alignment in alignments:
            if alignment is None:
                yield None
                continue

            # 1. Orientation & Finding Variants
            alignment, _ = GafAlignment.check_reverse(alignment, rgfa)
            
            logger.trace(f'Finding variants in {alignment.read_id}')
            # Use our new finding logic
            find_result = GAFReader.find_variants_in_alignment(
                alignment, variants, rgfa, variant_pointer
            )
            
            (variants_in_alignment, align_start_ref, new_pointer, 
             start_scaf_idx, end_scaf_idx, count_ext, count_sv) = find_result
            
            # Update global pointer for next read
            variant_pointer = new_pointer

            if not variants_in_alignment:
                yield (None, alignment, None)
                continue

            # 2. Build Reference Sequence & Calculate SV Attributes
            # We return the constructed sequence and the updated variants list
            reference_seq = self._calculate_sv_attributes(
                alignment, 
                variants_in_alignment, 
                rgfa, 
                start_scaf_idx
            )
            
            # 3. Interpolate External Variant Positions
            # SNPs don't need node summing; they just need anchor offsets
            self._interpolate_ext_positions(
                variants_in_alignment,
                alignment,
                rgfa,
                start_scaf_idx
            )
            
            yield (variants_in_alignment, alignment, reference_seq)


    def _alignments_to_reads(self, updated_variants):
        """Processes the identified variants on alignment and converts them in Read objects.

        Handles the conversion of variants detected into Read objects.
        Also performs the realignment of the reads with the alleles to get scores.

        Args:
            updated_variant: see the output of _update_variants_in_alignments

        Yields:
            Read: the Read object
        """

        cg_letter_to_op = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, 'X': 7, '=': 8}
        Alignment = namedtuple('Alignment', ['cigartuples', 'reference_start', 'query_sequence'])        # Class created to maintain compatibility with old code
        
        for result in updated_variants:

            # if no alignment was found 
            if result is None:
                yield None

            variants_in_alignment, alignment, reference = result

            # if no variants found in the alignment
            if variants_in_alignment is None:
                yield None

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
            )
            
            detected = self.detect_alleles_by_alignment(
                self._aligner,
                variants_in_alignment,
                0,                              # Here this has been hardcoded to 0. In original code, this was the index of the first variant (index in the big list of variants) in the read. But now we have new list of variants just for this alignment.
                processed_alignment,
                reference,
                self._overhang,
                self._is_custom_graph)
            for j, scores in detected:
                read.add_variant(variants_in_alignment[j].position_on_ref, scores)
            if read:  # At least one variant covered and detected
                yield read

    @staticmethod
    def realign(
        aligner: Realigner,
        variant: VcfVariant,
        read: Read,
        cigartuples: tuple,
        i: int,
        consumed: int,
        query_pos: int,
        reference: str,
        overhang: int,
        is_custom_graph: bool
    ):
        
        """Realigns a read to the reference and alternative alleles of a variant.

        Extracts the relevant query sequence and reference context (based on the 
        alignment path) to construct potential allele sequences. It then scores 
        each allele against the read using either simple edit distance (for 
        external variants) or a dedicated aligner (for SVs).

        Args:
            aligner: The aligner object used for scoring SVs.
            variant: The variant record to realign against.
            read: The AlignedSegment object (pysam).
            cigartuples: The cached `read.cigartuples` property (passed 
                explicitly to avoid expensive re-access).
            i: The index in `cigartuples` where the variant position occurs.
            consumed: The number of reference bases consumed by the CIGAR 
                operations up to index `i`.
            query_pos: The 0-based index in the query (read) sequence 
                corresponding to the variant's start position.
            reference: The specific sequence of the alignment path (constructed 
                from GAF nodes), NOT the entire chromosome sequence.
            overhang: The number of context bases to include on the left 
                and right of the variant for realignment.

        Returns:
            list[float]: A list of alignment scores, where the first element is 
            the Reference score, followed by scores for each Alternative allele.
            
            Returns `None` if the variant contains symbolic alleles (e.g., <DEL>).
        """
        # Do not process symbolic alleles like <DEL>, <DUP>, etc.
        if any([alt.startswith("<") for alt in variant.alternative_allele]):
            return None

        # There is a big difference between the previous implementation and what is needed.
        # In the previous code, the CIGAR is against the reference always and hence we need to realign only for the alternate alleles.
        # With GAF, the CIGAR is not always against the reference (sometimes it is not against ref or any of the alt and can be with a path that is not an allele traversal)
        # So we need to generalize the process to realign using the variant record and the cigar tuples.
        left_cigar, right_cigar = AlignmentReader.split_cigar(cigartuples, i, consumed)

        if not variant.is_sv():
            assert variant.state == 0
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
            for index, allele in enumerate([ref]+alts):
                scores.append(edit_distance(query, allele, aligner.bandwidth))
                
            return scores

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

        scores = []
        if is_custom_graph:
            # this means I can use the heuristics with the allele distance matrix
            
            closest_allele_idx = None   # this will store which allele is closest to the path the sequence aligned to.
            path_sequence = reference[variant.position : variant.position + variant.length_on_path] # this is the sequence of the path the variant aligned to.
            if path_sequence == ref_allele:
                closest_allele_idx = 0
            
            ref = left_overhang + ref_allele + right_overhang
            alts = []
            for idx, alt_allele in enumerate(variant.alternative_allele):
                if alt_allele != "*":
                    alt = left_overhang + alt_allele + right_overhang
                    if alt_allele == path_sequence:
                        closest_allele_idx = idx + 1
                else:
                    alt = left_overhang + right_overhang
                alts.append(alt)
            
            if closest_allele_idx != None:
                """
                if the allele with the exact seqeunce was not found,
                then it means that the path is a subsequence of the allele
                which can only happen with partial alignments
                """
                # this should not happen if all the alleles are single nodes!
                assert variant.state != 0
            
            '''
            New implementation:
            I have already calculate allele distance matrix which is the distance between allele i and allele j at this variant position. (Using edit distance)

            Now we use triangle inequality heuristics.
            d(A_i, R) <= |d(A_k, R) - d(A_k, A_i)| -> since d(A_k, R) <= b and d(A_k, A_i) <= 2b, this function is bounded by [0, 2b]

            We only want to do a realignment if d(A_i, R) < b.
            So we check if |d(A_k, R) - d(A_k, A_i)| < b and do realignment accordingly.
            '''
            if variant.state == 0:
                best_allele = ref if closest_allele_idx == 0 else alts[closest_allele_idx-1]
                best_score = aligner.get_distance(query, best_allele, 0)
                
                for idx, allele in enumerate([ref]+alts):
                    # TODO: There is the idea that if best_score == bandwidth,
                    # then that means that the best alignment is outside our alignment scope.
                    # So we can just automatically set all scores to bandwidth.
                    # PROBLEM: I don't trust aligners.
                    if idx == closest_allele_idx:
                        scores.append(best_score)
                        continue
                    idx1 = None
                    idx2 = None
                    if idx < closest_allele_idx:
                        idx1 = idx
                        idx2 = closest_allele_idx
                    else:
                        idx1 = closest_allele_idx
                        idx2 = idx
                    allele_to_allele_distance = variant.get_distance(idx1, idx2)
                    best_possible_score = abs(allele_to_allele_distance - best_score)
                    if best_possible_score >= aligner.bandwidth:
                        scores.append(aligner.bandwidth)
                    else:
                        scores.append(aligner.get_distance(query, allele, 0))
            else:
                # if its a partial alignment, then cannot apply the above heuristics
                for idx, allele in enumerate([ref]+alts):
                    scores.append(aligner.get_distance(query, allele, variant.state))
        else:
            # not a custom graph. 
            # so we will have multiple nodes for each alleles.
            # also the alignment path might not correspond to alleles.

            if variant.state == 0:
                for idx, allele in enumerate([ref]+alts):
                    if abs(len(query) - len(allele)) >= aligner.bandwidth:
                        # quick check of length difference.
                        scores.append(aligner.bandwidth)
                        continue
                    scores.append(aligner.get_distance(query, allele, variant.state))
            else:
                scores.append(aligner.get_distance(query, allele, variant.state))

        # Old implementation. Doing realignment for each allele.
        #for index, allele in enumerate([ref]+alts):
        #    if (abs(len(query) - len(allele)) > 5000 ) and (len(query)/len(allele) > 1.5 or len(query)/len(allele) < 1/1.5):
        #        # If the distance between the allele and query is too much, add a known high distance
        #        scores.append(1e8)
        #    else:
        #        scores.append(aligner.get_distance(query, allele))        
                
        return scores

    @staticmethod
    def detect_alleles_by_alignment(
        aligner: Realigner,
        variants: List[VcfVariant],
        j: int,
        read: Read,
        reference: str,
        overhang: int = 10,
        is_custom_graph: bool = False
    ):
        
        """Calculates the distance score between the read and the different alleles that are possible.

        Args:
            aligner: the aligner function - this can either be the edit_distance() function or a WaveFrontAligner object.
            variants: the variants present on the chromosome
            j: index of the first variant (in the variants list) to check
            read: the Read object of the alignment
            reference: the sequence of the refernce path (in the GAF file, it is the reference path) the read aligned to.
            mode:

        Returns:
            ReadSet: the ReadSet object containing the information about the variants found in each
                alignment and their distance scores from available alleles.
        """
        # Accessing read.cigartuples is expensive, do it only once
        cigartuples = read.cigartuples

        # For the same reason, the following check is here instad of
        # in the _usable_alignments method
        if not cigartuples:
            return
        for index, i, consumed, query_pos in _iterate_cigar(variants, j, read, cigartuples):
            scores = GAFReader.realign(
                aligner,
                variants[index],
                read,
                cigartuples,
                i,
                consumed,
                query_pos,
                reference,
                overhang,
                is_custom_graph
            )

            if scores is not None:
                yield (index, scores)


    def __enter__(self):
        return self

    def __exit__(self, *args):
        logger.debug("Closing GAFReader")
        self.close()

    def close(self):
        logger.debug("Deallocating WavefrontAligner")
        self._aligner.close()
        logger.debug("Closing GAFParser")
        self._reader.close()
