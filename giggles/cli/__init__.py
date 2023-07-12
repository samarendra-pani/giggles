import sys
import resource
import logging
from collections import defaultdict, namedtuple

from giggles.bam import (
    AlignmentFileNotIndexedError,
    EmptyAlignmentFileError,
    SampleNotFoundError,
    ReferenceNotFoundError,
)
from giggles.gaf import (
    AlignmentFileNotIndexedError,
    EmptyAlignmentFileError,
    ReferenceNotFoundError,
)
from giggles.variants import ReadSetReader, ReadSetError, GAFReader
from giggles.utils import IndexedFasta, FastaNotIndexedError, detect_file_format
from giggles.core import ReadSet
from giggles.vcf import VcfReader

logger = logging.getLogger(__name__)


class CommandLineError(Exception):
    """An anticipated command-line error occurred. This ends up as a user-visible error message"""


def alignment_reader(type, paths, reference, read_fasta, numeric_sample_ids, **kwargs):
    try:
        if type == "BAM":
            readset_reader = ReadSetReader(paths, reference, numeric_sample_ids, **kwargs)
        elif type == "GAF":
            readset_reader = GAFReader(paths, reference, read_fasta, numeric_sample_ids, **kwargs)
    except OSError as e:
        raise CommandLineError(e)
    except AlignmentFileNotIndexedError as e:
        raise CommandLineError(
            "The file '{}' is not indexed. Please create the appropriate BAM/CRAM "
            'index with "samtools index"'.format(e.args[0])
        )
    except EmptyAlignmentFileError as e:
        raise CommandLineError(
            "No reads could be retrieved from '{}'. If this is a CRAM file, possibly the "
            "reference could not be found. Try to use --reference=... or check your "
            "$REF_PATH/$REF_CACHE settings".format(e.args[0])
        )
    return readset_reader


class PhasedInputReader:
    def __init__(
        self,
        bam_or_gaf_paths,
        reference_fasta,
        gfa,
        read_fasta,
        numeric_sample_ids,
        **kwargs,  # passed to ReadSetReader constructor
    ):
        self._bam_paths, self._gaf_paths = self._split_input_file_list(bam_or_gaf_paths)
        if self._bam_paths and self._gaf_paths:
            raise CommandLineError("Unable to process both BAM and GAF files together. Please provide one.")
        if self._bam_paths:
            self._type="BAM"
            reference = reference_fasta
        else:
            self._type="GAF"
            reference = gfa
        logger.info("Detected %s file given as input..." %(self._type))
        self._numeric_sample_ids = numeric_sample_ids
        self._fasta = self._open_reference(reference_fasta) if reference_fasta else None

        self._readset_reader = alignment_reader(self._type, bam_or_gaf_paths, reference, read_fasta, numeric_sample_ids, **kwargs)
    
    def __enter__(self):
        return self

    def __exit__(self, *args):
        if self._fasta is not None:
            self._fasta.close()

    @property
    def has_vcfs(self):
        return bool(self._vcf_paths)

    @property
    def has_alignments(self) -> bool:
        """Whether any of the input files are BAM or CRAM"""
        return bool(self._bam_paths)

    @staticmethod
    def _split_input_file_list(paths):
        bams = []
        gafs = []
        for path in paths:
            try:
                file_format = detect_file_format(path)
            except OSError as e:
                raise CommandLineError(e)
            if file_format in ("BAM", "CRAM"):
                bams.append(path)
            elif file_format == "GAF":
                gafs.append(path)
            else:
                raise CommandLineError("Unable to determine type of input file {!r}".format(path))
        return bams, gafs

    @staticmethod
    def _open_reference(path):
        try:
            indexed_fasta = IndexedFasta(path)
        except OSError as e:
            raise CommandLineError("Error while opening FASTA reference file: {}".format(e))
        except FastaNotIndexedError as e:
            raise CommandLineError(
                "An index file (.fai) for the reference FASTA {!r} "
                "could not be found. Please create one with "
                '"samtools faidx".'.format(e)
            )
        return indexed_fasta

    def read(self, chromosome, variants, sample, haplotags, keep_untagged):
        """
        Return a pair (readset, vcf_source_ids) where readset is a sorted ReadSet.

        Set read_vcf to False to not read phased blocks from the VCFs
        """
        readset_reader = self._readset_reader
        for_sample = "for sample {!r} ".format(sample)
        logger.info("Reading alignments %sand detecting alleles ...", for_sample)
        reference = None
        if self._type == "BAM":
            try:
                reference = self._fasta[chromosome] if self._fasta else None
            except KeyError:
                raise CommandLineError(
                    "Chromosome {!r} present in VCF file, but not in the reference FASTA {!r}".format(
                        chromosome, self._fasta.filename
                    )
                )
            if reference == None:
                CommandLineError("No reference sequence found for Chromosomes {!r}. Please provide the reference file with the chromosomes.".format(chromosome))
        
        bam_sample = sample
        try:
            readset = readset_reader.read(chromosome, variants, bam_sample, reference)
        except SampleNotFoundError:
            logger.warning("Sample %r not found in any BAM/CRAM file.", bam_sample)
            readset = ReadSet()
        except ReadSetError as e:
            raise CommandLineError(e)
        except ReferenceNotFoundError:
            if chromosome.startswith("chr"):
                alternative = chromosome[3:]
            else:
                alternative = "chr" + chromosome
            message = "The chromosome {!r} was not found in the BAM/CRAM file.".format(chromosome)
            if readset_reader.has_reference(alternative):
                message += " Found {!r} instead".format(alternative)
            raise CommandLineError(message)

        new_readset = ReadSet()
        for read in readset:
            if not keep_untagged:
                if haplotags[read.name].hp != "none":
                    assert haplotags[read.name].hp == "H1" or haplotags[read.name].hp == "H2"
                    read.sort()
                    read.add_haplotag(haplotags[read.name].hp, haplotags[read.name].ps)
                    new_readset.add(read)
            else:
                read.sort()
                read.add_haplotag(haplotags[read.name].hp, haplotags[read.name].ps)
                new_readset.add(read)
        new_readset.sort()

        logger.info(
            "Found %d reads covering %d variants", len(new_readset), len(new_readset.get_positions())
        )
        return new_readset

def read_haplotags(file):
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

    logger.info("Reading Haplotag TSV File")
    Haplotag = namedtuple('Haplotag', ['hp', 'ps', 'chr'])
    out = defaultdict(lambda: Haplotag(hp="none", ps=-1, chr="none"))
    if file == None:
        return out
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
            out[rn] = Haplotag(hp=hp, ps=int(ps), chr=chr)
    return out



def log_memory_usage(include_children=False):
    if sys.platform == "linux":
        if include_children:
            memory_kb = (
                resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                + resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
            )
        else:
            memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum memory usage: %.3f GB", memory_kb / 1e6)
