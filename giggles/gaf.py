import os
from abc import ABC
from urllib.parse import urlparse
from typing import Iterator

import re
import pysam
import logging
import gzip
from collections import defaultdict, namedtuple
from dataclasses import dataclass
import pickle as pkl

logger = logging.getLogger(__name__)


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


def is_local(path):
    return urlparse(path).scheme == ""


def detect_gzip(path):
    with open(path, 'rb') as test_f:
        return (test_f.read(2) == b'\x1f\x8b')

class GafParser(ABC):
    pass

class GafAlignment:
    """
    Class to describe and work with GAF alignments
    """
    def __init__(self, line, source_id, fasta):
        self.read_id, self.q_len, self.q_start, self.q_end, self.orient, self.path, self.p_len, self.p_start, self.p_end, self.mapping_quality, self.tags, self.sequence = self.parseGafLine(line, fasta)
        if 'cg' in self.tags:
            self.cigar = self.tags['cg']
        else:
            self.cigar = self.tags['CG']
        self.path = list(filter(None, re.split('(>)|(<)', self.path)))
        self.source_id = source_id
        self.clip_start = False
        self.clip_end = False
        pass
    
    @staticmethod
    def parseGafLine(line, fasta):
        line = line.rstrip().split("\t")
        read_id = line[0]
        tags = {}
        for f in line[12:]:
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
        assert "tp" in tags, "No 'tp' tag to indicate primary alignment."
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


class SampleGafParser(GafParser):
    """
    Parsing the GAF file and extracting the alignment informartion.
    Since GAF files don't have sample specifications or multisample support, it will be assumed that all reads are for the same sample.
    """
    def __init__(
        self,
        path: str, *,
        reference: str,
        read_fasta: str = None,
        mapq = None,
        source_id: int = 0
    ):
        """
        path -- path to the GAF file
        reference -- rGFA for the realignment
        """
        reference = os.path.abspath(reference)
        self.source_id: int = source_id
        self._mapq = mapq
        path = os.path.abspath(path)
        logger.info("Reading GFA File.")
        self._reference = rGFA(reference_path=reference)
        if read_fasta != None:
            logger.info("Parsing FASTA File.")
            self._read_sequences = pysam.FastaFile(read_fasta)
        else:
            logger.info("FASTA file not given. Assuming the read sequences will be provided in GAF File.")
            self._read_sequences = None
        gzipped = detect_gzip(path)
        if gzipped:
            self._file = pysam.libcbgzf.BGZFile(path, "rb")
        else:
            self._file = open(path, "r")
        logger.info("Parsing GAF File.")
        self._index = self.process_index_file(path)

    def process_index_file(self, path):
        try:
            with open(path+".gai", 'rb') as f:
                return pkl.load(f)
        except FileNotFoundError:
            raise AlignmentFileNotIndexedError("No index file found for GAF file. Run gaftools scaffold-sort and create index.")

    def __call__(self, contig):
        self._contig_iter = contig
        assert self._reference.is_backbone(contig), "The contig specified is not a primary reference contig."
        return self

    def __iter__(self) -> Iterator[GafAlignment]:
        """
        Fetch GafAlignment from specified contig
        """
        offsets = self._index[self._contig_iter]    # Contains offset of first line of alignment and last line of alignment for a particular contig
        it = True
        file = self._file
        file.seek(offsets[0])
        while it == True:
            line = file.readline()
            if not line:
                break
            a = GafAlignment(line, self.source_id, self._read_sequences)
            if a.mapping_quality < self._mapq:
                continue
            if a.tags['tp'] != "P":
                continue
            if a.tags['sn'] != self._contig_iter:
                assert a.tags['sn'] == 'unknown', "GAF is not properly sorted."
                continue
            if a.tags['iv'] == 1:
                continue
            # TODO: Multiple alignments with same read ids have to be checked later and one selected.
            if file.tell() == offsets[1]:
                it = False
            
            yield a



class rGFA:
    """
    Parsing the reference GFA file to find the main reference backbone.
    Also store information about the reference sequence.
    """
    def __init__(self, reference_path) -> None:
        
        
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
        Node = namedtuple("Node", ['sequence', 'start', 'contig', 'tags'])
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
            node_dict[node_id] = Node(node_seq, node_start, node_contig, tags)
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