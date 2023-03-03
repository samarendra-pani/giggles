import sys
import resource
import logging

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


def open_readset_reader(*args, **kwargs):
    try:
        readset_reader = ReadSetReader(*args, **kwargs)
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

def open_gaf_reader(*args, **kwargs):
    try:
        readset_reader = GAFReader(*args, **kwargs)
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
        else:
            self._type="GAF"
        logger.info("Detected %s file given as input..." %(self._type))
        self._numeric_sample_ids = numeric_sample_ids
        self._fasta = self._open_reference(reference_fasta) if reference_fasta else None

        if self._type == "BAM":
            self._readset_reader = open_readset_reader(self._bam_paths, reference_fasta, numeric_sample_ids, **kwargs)
        else:
            self._readset_reader = open_gaf_reader(self._gaf_paths, gfa, read_fasta, numeric_sample_ids, **kwargs)
    
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

    def read(self, chromosome, variants, sample, *, regions=None):
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
        elif self._type == "GAF":
            reference = self._readset_reader._reader._reference.get_sequence(chromosome)

        if reference == None:
            CommandLineError("No reference sequence found for Chromosomes {!r}. Please provide the reference file with the chromosomes.".format(chromosome))
        bam_sample = sample
        try:
            readset = readset_reader.read(chromosome, variants, bam_sample, reference, regions)
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

        for read in readset:
            read.sort()
        readset.sort()

        logger.info(
            "Found %d reads covering %d variants", len(readset), len(readset.get_positions())
        )
        return readset

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
