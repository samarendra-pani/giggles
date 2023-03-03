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
    def __init__(self, line, offset, source_id):
        self.read_id, self.q_len, self.q_start, self.q_end, self.orient, self.path, self.p_len, self.p_start, self.p_end, self.mapping_quality = self.parseGafLine(line)
        self.path = list(filter(None, re.split('(>)|(<)', self.path)))
        self.offset = offset
        self.source_id = source_id
        self.sequence = None
        self.tags = None
        self.clip_start = False
        self.clip_end = False
        pass
    
    @staticmethod
    def parseGafLine(line):
        line = line.rstrip().split("\t")
        return line[0], int(line[1]), int(line[2]), int(line[3]), line[4], line[5], int(line[6]), int(line[7]), int(line[8]), int(line[11])

    def detect_chromosome(self, gfa_ref):
        """
        Take a dictionary as input. The dictionary is SampleGafParser._reference._nodes
        """
        primary_contigs = gfa_ref._primary_contigs
        path = self.path
        chr = None
        for n, i in enumerate(path):
            if i == "<" or i == ">":
                continue
            found = False
            for contig in primary_contigs:
                if i in gfa_ref._nodes[contig][0]:
                    if chr != None:
                        assert contig == chr, "Read has been mapped to two reference contigs."
                    chr = contig
                    found = True
                    break
            if (n == 2) and not found:
                self.clip_start = True
            if (n == len(path)) and not found:    
                self.clip_end = True
        
        return chr

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
            self._file = pysam.libcbgzf.BGZFile(path, "rt")
        else:
            self._file = open(path, "r")
        logger.info("Parsing GAF File.")
        self._alignments = self.parseGAF(mapq)
        
    def parseGAF(self, mapq):
        
        file = self._file
        gfa_ref = self._reference
        alignments = {}
        while True:
            offset = file.tell()
            line = file.readline()
            if not line:
                break
            fields = line.split("\t")
            read_id = fields[0]
            a = GafAlignment(line, offset, self.source_id)
            if a.mapping_quality < mapq:
                continue
            chr = a.detect_chromosome(gfa_ref)
            # This filters out reads not mapped to any primary reference contigs
            if chr == None:
                continue
            tp = None
            for f in fields[12:]:
                if not f.startswith("tp:A:"):
                    continue
                tp = f[5]
            if tp != "P":
                continue
            if chr not in alignments:
                alignments[chr] = {}
                alignments[chr][read_id] = a
                continue
            # This compares alignments with same read ids and selects one based on mapping quality or number of nodes.
            if read_id in alignments[chr]:
                keep = a.compare(alignments[chr][read_id])
                if keep:
                    alignments[chr][read_id] = a
                continue
            alignments[chr][read_id] = a
            
        return alignments

    def __call__(self, contig):
        self._contig_iter = contig
        assert self._reference.is_primary(contig), "The contig specified is not a primary reference contig."
        return self

    def __iter__(self) -> Iterator[GafAlignment]:
        """
        Fetch GafAlignment from specified contig
        """
        for read_id, alignment in self._alignments[self._contig_iter].items():
            assert read_id == alignment.read_id
            off = alignment.offset
            self._file.seek(off)
            line = self._file.readline()
            tags = {}
            seq = None
            for f in line.split("\t")[12:]:
                f = f.split(":")
                assert len(f) ==  3, "The tag provided in read %s is not in correct format."%(read_id)
                if f[1] == "i":
                    f[2] = int(f[2])
                elif f[1] == "f":
                    f[2] = float(f[2])
                tags[f[0]] = f[2]
            assert "cg" in tags or "CG" in tags, "No CIGAR string in read %s. Provide the CIGAR string with 'CG' or 'cg' tag."%(read_id)
            # The sequence will be searched in the RS tag or rs tag
            if self._read_sequences == None:
                assert "rs" in tags or "RS" in tags, "No Read Sequence in read %s. Provide the CIGAR string with 'RS' or 'rs' tag."%(read_id)
                if "rs" in tags:
                    seq = tags.pop("rs")
                else:
                    seq = tags.pop("RS")
            else:
                seq = self._read_sequences.fetch(read_id)
            alignment.set_tags(tags)
            alignment.set_sequence(seq)
            yield alignment







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
        self._nodes = self.parse_gfa_file(file)
        file.close()
        self._contigs = list(self._nodes.keys())
        self._primary_contigs = [x for x in self._contigs if self._nodes[x][1] == 0]


    @staticmethod
    def parse_gfa_file(file):
        Node = namedtuple("Node", ['sequence', 'start'])
        node_dict = {}
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
            node_contig = tags["SN"]
            node_start = tags["SO"]
            node_rank = tags["SR"]
            try:
                node_dict[node_contig][0][node_id] = Node(node_seq, node_start)
                assert (node_rank == node_dict[node_contig][1])
            except:
                node_dict[node_contig] = [{}, node_rank]
                node_dict[node_contig][0][node_id] = Node(node_seq, node_start)
        
        #Sorting the nodes.
        #TODO: Put a check if sorting is required
        for contig in node_dict.keys():
            node_dict[contig][0] = dict(sorted(node_dict[contig][0].items(), key=lambda x: x[1].start))
            
        return node_dict
    
    def is_primary(self, contig):
        assert (contig in self._contigs)
        return self._nodes[contig][1] == 0
    
    def get_sequence(self, contig):
        if self._nodes[contig][1] != 0:
            error = "The chromosome specified is not in the reference. The reference contigs are: "
            for i in self._nodes.keys():
                error += "\n%s"%(str(i))
            raise Exception(error)
        nodes = self._nodes[contig][0]
        sequence = ""
        pos = 0
        for _, value in nodes.items():
            if value.start != pos:
                raise Exception('The reference nodes either have overlap or have empty space. Cannot reconstruct reference sequence for chromosome %s'%(contig))
            sequence += value.sequence
            pos += len(value.sequence)
        
        return sequence