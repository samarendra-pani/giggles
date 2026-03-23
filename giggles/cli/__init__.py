# Code modified from WhatsHap (https://github.com/whatshap/whatshap)

import sys
import resource
from giggles.logger import logger

from giggles.variants import GAFReader
from giggles.core import ReadSet


class ReadSetCreator:
    def __init__(
        self,
        alignment_files,
        gfa,
        read_fasta_files,
        **kwargs,  # passed to GAFReader constructor
    ):
        self.readset_reader = GAFReader(alignment_files, gfa, read_fasta_files, **kwargs)

    def __exit__(self, *args):
        self.readset_reader.close()

    def __enter__(self):
        return self

    def read(self, chromosome, variants, haplotags, keep_untagged):
        """
        Returns a ReadSet object with haplotags
        """

        logger.info("Reading alignments and detecting alleles.")
        readset = self.readset_reader.read(chromosome, variants)
        if readset is None:
            readset = ReadSet()
        
        new_readset = ReadSet()
        for read in readset:
            read_id = (read.source_id, read.name)
            try:
                haplotag = haplotags[read_id]
                if not keep_untagged:
                    if haplotag.hp != 'none':
                        assert haplotag.hp == 'H1' or haplotag.hp == 'H2'
                        read.sort()
                        if haplotag.hp == 'H1':
                            read.add_haplotag(False)
                        else:
                            read.add_haplotag(True)
                        read.add_phaseset(haplotag.ps)
                        new_readset.add(read)
                else:
                    read.sort()
                    read.add_haplotag(haplotag.hp)
                    read.add_phaseset(haplotag.ps)
                    new_readset.add(read)
            except KeyError:
                logger.warning(f'Could not find haplotag for read {read_id[1]} from source file {read_id[0]}.')
                if keep_untagged:
                    read.sort()
                    new_readset.add(read)
            except TypeError:
                assert haplotags is None
                assert keep_untagged
                read.sort()
                new_readset.add(read)

        new_readset.sort()

        return new_readset


class Haplotag:
    def __init__(self, hp, ps, chr):
        self.hp = hp
        self.ps = ps
        self.chr = chr


def read_haplotags(files):
    """
    Function to read the haplotag file.
    The file should be tab-separated with the following column information:
    Column 1 - readname
    Column 2 - haplotype (H1 or H2 or none)
    Column 3 - phaseset (check whatshap documentation for information on PS)
    Column 4 - chromosome
    Assuming that any title or comment lines start with '#'.

    This is the standard output for `whatshap haplotag`. Check documentation for more information.
    """

    haplotags = {}
    if files == None:
        logger.info("No haplotag files provided.")
        return None
    for index, file in enumerate(files):
        with open(file, 'r') as f:
            while True:
                line = f.readline()
                if not line:
                    break
                if line[0] == "#":
                    continue
                rn, hp, ps, chr = line.rstrip().split('\t')[0:4]
                if ps == 'none':
                    ps = -1
                assert hp in ['H1', 'H2', 'none']
                haplotags[(index, rn)] = Haplotag(hp=hp, ps=int(ps), chr=chr)
    logger.info(f"Found {len(haplotags)} haplotags.")
    return haplotags



def log_memory_usage(include_children=False):
    if sys.platform == "linux":
        if include_children:
            memory_kb = (
                resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                + resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
            )
        else:
            memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info(f"Maximum memory usage: {memory_kb / 1e6:.3f} GB")
