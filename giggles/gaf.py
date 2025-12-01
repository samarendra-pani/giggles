import os
from abc import ABC
from urllib.parse import urlparse
from typing import Iterator
from .utils import reverse_complement, reverse_cigar
from .logger import logger, warn_once


import re
import pysam
import gzip
from dataclasses import dataclass
import pickle as pkl


@dataclass
class AlignmentWithSourceID:
    source_id: int
    bam_alignment: pysam.AlignedSegment


def is_local(path):
    return urlparse(path).scheme == ""


def detect_gzip(path):
    with open(path, 'rb') as test_f:
        return (test_f.read(2) == b'\x1f\x8b')


'''
Future Idea:
    - Work with GAF alignments with a reduced representation.
    - Use the BO and NO index and collapse completely encompassed bubbles into single nodes.
    - Store the bare minimum values for downstream processing.
'''
class GafAlignment:
    """
    Class to describe and work with GAF alignments
    """
    def __init__(self, line, source_id, fasta):
        self.read_id, self.q_len, self.q_start, self.q_end, self.orient, self.path, self.p_len, self.p_start, self.p_end, self.mapping_quality, self.tags, self.sequence = self.parseGafLine(line, fasta)
        if 'cg' in self.tags:
            self.cigar = self.tags.pop('cg')
        else:
            self.cigar = self.tags.pop('CG')
        self.path = list(filter(None, re.split('(>)|(<)', self.path)))
        self.source_id = source_id
        self.clip_start = False
        self.clip_end = False
        
    # Checking if the alignment is in reverse direction.
    @staticmethod
    def check_reverse(alignment, rgfa):
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
            alignment.sequence = reverse_complement(alignment.sequence)
            alignment.path = new_alignment
            alignment.cigar = reverse_cigar(alignment.cigar)
        
        return alignment, reverse

    @staticmethod
    def parseGafLine(line, fasta):
        try:
            line = line.rstrip().split("\t")
        except TypeError:
            line = line.decode('utf8').rstrip().split("\t")
        read_id = line[0]
        tags = {}
        for f in line[12:]:
            if '::' in f:
                # Right now making an exception for the 'ds' tag
                assert f.startswith("ds")
                f = f.split("::")
                tag_name = f[0].split(':')[0]
                tag_type = f[0].split(':')[1]
                if tag_type == "i":
                    f[1] = int(f[1])
                elif tag_type == "f":
                    f[1] = float(f[1])
                tags[tag_name] = f[1]
                continue 
            f = f.split(":")
            assert len(f) ==  3, "The tag provided in read %s is not in correct format."%(read_id)
            if f[1] == "i":
                f[2] = int(f[2])
            elif f[1] == "f":
                f[2] = float(f[2])
            tags[f[0]] = f[2]
        assert "cg" in tags or "CG" in tags, "No CIGAR string in read %s. Provide the CIGAR string with 'CG' or 'cg' tag."%(read_id)
        assert "sn" in tags, "No 'sn' tag to indicate contig it is aligned to. Use gaftools scaffold-sort to sort it and automatically add tag."
        assert "iv" in tags, "No 'iv' tag to indicate presence of inversion. Use gaftools scaffold-sort to sort it and automatically add tag."
        if "tp" not in tags:
            warn_once(logger, "No 'tp' tag to indicate primary alignment. Assuming the alignment is primary.")
            tags["tp"] = "P"
        # The sequence will be searched in the RS tag or rs tag
        if fasta == None:
            assert "rs" in tags or "RS" in tags, "No Read Sequence in read %s. Provide the CIGAR string with 'RS' or 'rs' tag."%(read_id)
            if "rs" in tags:
                seq = tags.pop("rs")
            else:
                seq = tags.pop("RS")
        else:
            seq = fasta.fetch(region=read_id)
        return line[0], int(line[1]), int(line[2]), int(line[3]), line[4], line[5], int(line[6]), int(line[7]), int(line[8]), int(line[11]), tags, seq

    @staticmethod
    def get_alignment_start_on_ref(alignment, rgfa):
        """Calculates the start position of the alignment on the reference.

        Iterates through the alignment path to find the first reference node
        and calculates the coordinate based on the distance from the path start.

        Args:
            alignment (GafAlignment): The alignment object
            rgfa: The rGFA graph object.

        Returns:
            tuple[int, int, int]: A tuple containing:
                - Start position on reference (0-based, inclusive).
                - Index of the start reference node in the path.
                - Index of the start scaffold node in the path.

            Returns (None, None, None) if no reference node is found.
        """
        start_on_path = alignment.p_start
        path = alignment.path
        count = -1
        start_node_on_ref = None
        found_scaffold = False
        found_ref = False
        scaffold_count = None
        ref_count = None
        for node in path:
            count += 1
            if node in ['<', '>']:
                continue
            if not found_scaffold:
                if rgfa.get_node(node).tags['NO'] == 0:
                    found_scaffold = True
                    scaffold_count = count
            if not found_ref:
                if rgfa.get_node(node).rank == 0:
                    found_ref = True
                    start_node_on_ref = node
                    ref_count = count
            if found_scaffold and found_ref:
                break
        if start_node_on_ref == None:
            return None, None, None
        assert ref_count > 0
        if ref_count == 1:
            # first node in the path is a reference node
            return rgfa.get_node(start_node_on_ref).start + start_on_path, ref_count, scaffold_count
        else:
            # first node in the path is not a reference node
            return rgfa.get_node(start_node_on_ref).start, ref_count, scaffold_count
    
    @staticmethod
    def get_alignment_end_on_ref(alignment, rgfa):
        """Calculates the end position of the alignment on the reference.

        Iterates backwards through the alignment path to find the last reference 
        node and calculates the coordinate based on the distance from the path end.

        Args:
            alignment (GafAlignment): The alignment object
            rgfa: The rGFA graph object.

        Returns:
            tuple[int, int, int]: A tuple containing:
                - End position on reference (0-based, exclusive).
                - Index of the end reference node in the path.
                - Index of the end scaffold node in the path.

            Returns (None, None, None) if no reference node is found.
        """
        distance_from_end = alignment.p_len - alignment.p_end
        path = alignment.path
        count = -1
        end_node_on_ref = None
        found_scaffold = False
        found_ref = False
        scaffold_count = None
        ref_count = None
        for node in path[::-1]:
            count += 1
            if node in ['<', '>']:
                continue
            if not found_scaffold:
                if rgfa.get_node(node).tags['NO'] == 0:
                    found_scaffold = True
                    scaffold_count = count
            if not found_ref:
                if rgfa.get_node(node).rank == 0:
                    found_ref = True
                    end_node_on_ref = node
                    ref_count = count
            if found_ref and found_scaffold:
                break
        if end_node_on_ref == None:
            return None, None, None
        if ref_count != None:
            ref_count = len(path) - ref_count - 1
        if scaffold_count != None:
            scaffold_count = len(path) - scaffold_count - 1
        if ref_count == len(path) - 1:
            # last node in the path is a reference node
            return rgfa.get_node(end_node_on_ref).start + len(rgfa.get_node(end_node_on_ref).sequence) - distance_from_end, ref_count, scaffold_count
        else:
            # last node in the path is not a reference node
            return rgfa.get_node(end_node_on_ref).start + len(rgfa.get_node(end_node_on_ref).sequence), ref_count, scaffold_count
    

    def compare(self, alignment):
        """
        Compare this alignment to another alignment object
        """
        # First check number of nodes covered by them
        if len(self.path) > len(alignment.path):
            return True
        elif len(self.path) < len(alignment.path):
            return False
        
        # Check length of alignment
        if self.q_len > alignment.q_len:
            return True
        elif self.q_len < alignment.q_len:
            return False

        if self.mapping_quality > alignment.mapping_quality:
            return True
        else:
            return False

    def __repr__(self):
        return f"GafAlignment(read_id={self.read_id}, q_len={self.q_len}, q_start={self.q_start}, q_end={self.q_end}, orient={self.orient}, path={''.join(self.path)}, p_len={self.p_len}, p_start={self.p_start}, p_end={self.p_end}, tags={self.tags}, cigar='{self.cigar}')"

    def set_tags(self, tags):
        self.tags = tags

    def has_tag(self, tag):
        return tag in self.tags

    def get_tag(self, tag):
        return self.tags[tag]

    def set_sequence(self, sequence):
        self.sequence = sequence

    def __del__(self):
        pass


class GafParser:
    """
    Parsing the GAF file and extracting the alignment informartion.
    Since GAF files don't have sample specifications or multisample support, it will be assumed that all reads are for the same sample.
    """
    def __init__(
        self,
        alignment_files,
        reference,
        read_fasta_files,
        mapq,
    ):
        """
        path -- path to the GAF file
        reference -- rGFA for the realignment
        """
        self._mapq = mapq
        self._files = []
        self._fastas = []
        self._indexes = []
        self._alignment_files = alignment_files
        alignment_files = [os.path.abspath(f) for f in alignment_files]
        read_fasta_files = [os.path.abspath(f) for f in read_fasta_files]
        self._reference = reference
        assert len(alignment_files) == len(read_fasta_files)
        for path in alignment_files:
            if detect_gzip(path):
                self._files.append(pysam.libcbgzf.BGZFile(path, "rb"))
            else:
                self._files.append(open(path, "r"))
            self._indexes.append(self.process_index_file(path))
        for path in read_fasta_files:
            self._fastas.append(pysam.FastaFile(path))
        logger.info("Completed initializing GAF files, along with their indexes and read FASTA files.")
 
    def process_index_file(self, path):
        try:
            logger.info(f"Processing index file for {path}")
            with open(path+".gsi", 'rb') as f:
                return pkl.load(f)
        except FileNotFoundError:
            raise Exception(f"No index file found for {path}. Run gaftools sort and create index.")

    def __call__(self, contig):
        self._contig_iter = contig
        assert self._reference.is_backbone(contig), "The contig specified is not a primary reference contig."
        return self

    def __iter__(self) -> Iterator[GafAlignment]:
        """
        Fetch GafAlignment from specified contig
        """
        for source_id, alignment_file in enumerate(self._files):
            try:
                # Contains offset of first line of alignment and last line of alignment for a particular contig
                offsets = self._indexes[source_id][self._contig_iter]
            except KeyError:
                logger.debug(f"No alignments found for contig {self._contig_iter} in file {self._alignment_files[source_id]}.")
                yield None
                continue
            alignment_file.seek(offsets[0])
            iterator = True
            while iterator:
                line = alignment_file.readline()
                if not line:
                    break
                if alignment_file.tell() == offsets[1]:
                    iterator = False
                a = GafAlignment(line, source_id, self._fastas[source_id])
                if a.mapping_quality < self._mapq:
                    continue
                if a.tags['tp'] != "P":
                    continue
                if a.tags['sn'] != self._contig_iter:
                    assert a.tags['sn'] == 'unknown', "GAF is not properly sorted."
                    continue
                # TODO: What to do with this inversion case?
                if a.tags['iv'] == 1:
                    continue
                yield a

    def __exit__(self):
        self.close()

    def close(self):
        for f in self._files:
            f.close()
        for f in self._fastas:
            f.close()
        for f in self._indexes:
            del f

class Node:
    def __init__(self, sequence, start, contig, rank, tags):
        self.sequence = sequence
        self.start = start
        self.contig = contig
        self.rank = rank
        self.tags = tags

    def __repr__(self):
        return "Node(sequence={}, start={}, contig={}, rank = {}, tags={})".format(self.sequence, self.start, self.contig, self.rank, self.tags)

class rGFA:
    """
    Parsing the reference GFA file to find the main reference backbone.
    Also store information about the reference sequence.
    """
    def __init__(self, reference_path) -> None:
        
        logger.info('Reading rGFA file.')
        gzipped = None
        with open(reference_path, 'rb') as test_f:
            gzipped = (test_f.read(2) == b'\x1f\x8b')
        if gzipped:
            file = gzip.open(reference_path, "rt")
        else:
            file = open(reference_path, "r")
        self.parse_gfa_file(file)
        file.close()
        

    def parse_gfa_file(self, file):
        node_dict = {}
        ref_contig_nodes = {}   # Node IDs (sequential) for reference backbone contigs 
        for line in file:
            if line[0] != "S":
                continue
            line = line.rstrip().split("\t")
            node_id = line[1]
            node_seq = line[2]
            tags = {}
            for i in line[3:]:
                i = i.split(":")
                if len(i) != 3:
                    continue
                if i[1] == "i":
                    tags[i[0]] = int(i[2])
                else:
                    tags[i[0]] = i[2]
            # Assuming that these three tags are present
            node_contig = tags.pop("SN")
            node_start = tags.pop("SO")
            node_rank = tags.pop("SR")
            node_dict[node_id] = Node(node_seq, node_start, node_contig, node_rank, tags)
            if node_rank == 0:
                try:
                    ref_contig_nodes[node_contig].append(node_id)
                except KeyError:
                    ref_contig_nodes[node_contig] = [node_id]

        #Sorting the nodes.
        #TODO: Put a check if sorting is required
        for contig in ref_contig_nodes.keys():
            ref_contig_nodes[contig] = sorted(ref_contig_nodes[contig], key=lambda x: node_dict[x].start)
             
        self._nodes = node_dict
        self._ref_contig_nodes = ref_contig_nodes

    def is_backbone(self, contig):
        return contig in self._ref_contig_nodes.keys()
    
    def get_backbone_sequence(self, contig):
        if contig not in self._ref_contig_nodes.keys():
            error = "The chromosome specified is not in the reference. The reference contigs are: "
            for i in self._ref_contig_nodes.keys():
                error += "\n%s"%(str(i))
            raise Exception(error)
        nodes = self._ref_contig_nodes[contig]
        sequence = ""
        pos = 0
        for node in nodes:
            node = self._nodes[node]
            if node.start != pos:
                raise Exception('The reference nodes either have overlap or have empty space. Cannot reconstruct reference sequence for chromosome %s'%(contig))
            sequence += node.sequence
            pos += len(node.sequence)
        
        return sequence
    
    def get_node(self, node_id):
        return self._nodes[node_id]