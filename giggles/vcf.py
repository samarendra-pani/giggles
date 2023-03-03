"""
Functions for reading VCFs.
"""
import os
import re
import sys
import math
import logging
import itertools
from dataclasses import dataclass
from abc import ABC, abstractmethod
from os import PathLike
from typing import List, Sequence, Dict, Tuple, Iterable, Optional, Union, TextIO, Iterator

from pysam import VariantFile, VariantHeader, VariantRecord

from .core import (
    PhredGenotypeLikelihoods,
    Genotype,
    binomial_coefficient,
)
from .utils import warn_once

logger = logging.getLogger(__name__)


class VcfError(Exception):
    pass


class VcfNotSortedError(VcfError):
    pass


class PloidyError(VcfError):
    pass


class VcfIndexMissing(VcfError):
    pass


class VcfInvalidChromosome(VcfError):
    pass


@dataclass
class VariantCallPhase:
    block_id: int  # numeric id of the phased block
    phase: Tuple[int, ...]  # alleles representing the phasing. (1, 0) is 1|0
    quality: Optional[int]


class VcfVariant:
    """A variant in a VCF file (not to be confused with core.Variant)"""

    __slots__ = ("position", "reference_allele", "alternative_allele", "allele_origin", "allele_traversal", "start_node", "end_node")

    def __init__(self, position: int, reference_allele: str, alternative_allele: tuple, allele_origin: list, allele_traversal: tuple):
        """
        Multi-ALT sites are not modelled.
        """
        self.position = position
        self.reference_allele = reference_allele
        self.alternative_allele = alternative_allele
        self.allele_origin = allele_origin
        self.allele_traversal = allele_traversal
        self.start_node = None
        self.end_node = None
        self.extract_edge_nodes()
    
    def extract_edge_nodes(self):
        for path in self.allele_traversal:
            p = list(filter(None, re.split('(>)|(<)', path)))
            if p[0] == ">" and p[-2] == ">":
                self.start_node = p[1]
                self.end_node = p[-1]
            else:
                raise Exception("Found Bubble with incorrect orientaion for the start and end node.")

    def __repr__(self):
        return "VcfVariant({}, {!r}, {!r}, {!r})".format(
            self.position, self.reference_allele, self.alternative_allele, self.allele_origin
        )

    def __hash__(self):
        return hash((self.position, self.reference_allele, self.alternative_allele))

    def __eq__(self, other):
        return (
            (self.position == other.position)
            and (self.reference_allele == other.reference_allele)
            and (self.alternative_allele == other.alternative_allele)
        )

    def __lt__(self, other):
        return (self.position, self.reference_allele, self.alternative_allele) < (
            other.position,
            other.reference_allele,
            other.alternative_allele,
        )

    def is_snv(self, ix) -> bool:
        return (self.reference_allele != self.alternative_allele[ix]) and (
            len(self.reference_allele) == len(self.alternative_allele[ix]) == 1
        )
       


class GenotypeLikelihoods:
    __slots__ = "log_prob_genotypes"

    def __init__(self, log_prob_genotypes: List[float]):
        """Likelihoods of all genotypes to be given as log10 of
        the original probability."""
        self.log_prob_genotypes = log_prob_genotypes

    def __repr__(self):
        return "GenotypeLikelihoods({})".format(self.log_prob_genotypes)

    def __eq__(self, other):
        if other is None:
            return False
        if self.log_prob_genotypes is None and other.log_prob_genotypes is None:
            return True
        return self.log_prob_genotypes == other.log_prob_genotypes

    def log10_probs(self) -> List[float]:
        return self.log_prob_genotypes

    def log10_prob_of(self, genotype_index: int) -> float:
        return self.log10_probs()[genotype_index]

    def as_phred(self, ploidy: int = 2, regularizer: float = None) -> PhredGenotypeLikelihoods:
        if regularizer is None:
            # shift log likelihoods such that the largest one is zero
            m = max(self.log_prob_genotypes)
            return PhredGenotypeLikelihoods(
                [round((prob - m) * -10) for prob in self.log_prob_genotypes], ploidy=ploidy
            )
        else:
            p = [10 ** x for x in self.log_prob_genotypes]
            s = sum(p)
            p = [x / s + regularizer for x in p]
            m = max(p)
            return PhredGenotypeLikelihoods(
                [round(-10 * math.log10(x / m)) for x in p], ploidy=ploidy
            )


class VariantTable:
    """
    For a single chromosome, store variants and their genotypes.
    Each row of this table contains a variant, each column
    contains the genotypes of a single sample.

    chromosome -- chromosome name
    samples -- list of sample names
    """

    def __init__(self, chromosome: str, samples: List[str]):
        self.chromosome = chromosome
        self.samples = samples
        self.genotypes: List[List[Genotype]] = [[] for _ in samples]
        self.phases: List[List[Optional[VariantCallPhase]]] = [[] for _ in samples]
        self.genotype_likelihoods: List[List[Optional[GenotypeLikelihoods]]] = [[] for _ in samples]
        self.variants: List[VcfVariant] = []
        self._sample_to_index = {sample: index for index, sample in enumerate(samples)}

    def __len__(self) -> int:
        return len(self.variants)

    def add_variant(
        self,
        variant: VcfVariant,
        genotypes: Sequence[Genotype],
        phases: Sequence[Optional[VariantCallPhase]],
        genotype_likelihoods: Sequence[Optional[GenotypeLikelihoods]],
    ) -> None:
        """Add a row to the table"""
        if len(genotypes) != len(self.genotypes):
            raise ValueError("Expecting as many genotypes as there are samples")
        if len(phases) != len(self.phases):
            raise ValueError("Expecting as many phases as there are samples")
        self.variants.append(variant)
        for i, genotype in enumerate(genotypes):
            assert isinstance(genotype, Genotype)
            self.genotypes[i].append(genotype)
        for i, phase in enumerate(phases):
            self.phases[i].append(phase)
        for i, gl in enumerate(genotype_likelihoods):
            self.genotype_likelihoods[i].append(gl)

    def genotypes_of(self, sample: str) -> List[Genotype]:
        """Retrieve genotypes by sample name"""
        return self.genotypes[self._sample_to_index[sample]]

    def set_genotypes_of(self, sample: str, genotypes: List[Genotype]) -> None:
        """Set genotypes by sample name"""
        assert len(genotypes) == len(self.variants)
        self.genotypes[self._sample_to_index[sample]] = genotypes

    def genotype_likelihoods_of(self, sample: str) -> List[Optional[GenotypeLikelihoods]]:
        """Retrieve genotype likelihoods by sample name"""
        return self.genotype_likelihoods[self._sample_to_index[sample]]

    def set_genotype_likelihoods_of(
        self, sample: str, genotype_likelihoods: List[Optional[GenotypeLikelihoods]]
    ) -> None:
        """Set genotype likelihoods by sample name"""
        assert len(genotype_likelihoods) == len(self.variants)
        self.genotype_likelihoods[self._sample_to_index[sample]] = genotype_likelihoods

    def phases_of(self, sample: str) -> List[Optional[VariantCallPhase]]:
        """Retrieve phases by sample name"""
        return self.phases[self._sample_to_index[sample]]

    def num_of_blocks_of(self, sample: str) -> int:
        """ Retrieve the number of blocks of the sample"""
        return len(
            set(i.block_id for i in self.phases[self._sample_to_index[sample]] if i is not None)
        )

    def id_of(self, sample: str) -> int:
        """Return a unique int id of a sample given by name"""
        return self._sample_to_index[sample]


class MixedPhasingError(Exception):
    pass


class VcfReader:
    """
    Read a VCF file chromosome by chromosome.
    """

    def __init__(
        self,
        path: Union[str, PathLike],
        bam_samples: List[str] = None,
        indels: bool = False,
        phases: bool = False,
        genotype_likelihoods: bool = False,
        ignore_genotypes: bool = False,
        ploidy: int = None,
    ):
        """
        path -- Path to VCF file
        indels -- Whether to include also insertions and deletions in the list of
            variants.
        ignore_genotypes -- In case of genotyping algorithm, no genotypes may be given in
                                vcf, so ignore all genotypes
        ploidy -- Ploidy of the samples
        """
        # TODO Always include deletions since they can 'overlap' other variants
        self._indels = indels
        self._vcf_reader = VariantFile(os.fspath(path))
        self._path = path
        self._phases = phases
        self._genotype_likelihoods = genotype_likelihoods
        self._ignore_genotypes = ignore_genotypes
        self.samples = bam_samples  # intentionally public   # Make them into BAM File samples
        self.ploidy = ploidy
        logger.debug("Found %d sample(s) in the BAM file.", len(self.samples))      # BAM File samples

    def __enter__(self):
        return self

    def __exit__(self, *args):
        # follows same structure as for ReadSetReader
        self.close()

    def close(self):
        self._vcf_reader.close()

    @property
    def path(self) -> str:
        return self._vcf_reader.filename.decode()

    def _fetch(self, chromosome: str, start: int = 0, end: Optional[int] = None):
        try:
            records = self._vcf_reader.fetch(chromosome, start=start, stop=end)
        except ValueError as e:
            if "invalid contig" in e.args[0]:
                raise VcfInvalidChromosome(e.args[0]) from None
            elif "fetch requires an index" in e.args[0]:
                raise VcfIndexMissing(
                    "{} is missing an index (.tbi or .csi)".format(self._path)
                ) from None
            else:
                raise
        return records

    def fetch(self, chromosome: str, start: int = 0, end: Optional[int] = None) -> VariantTable:
        """
        Fetch records from a single chromosome, optionally restricted to a single region.

        Return a VariantTable object.
        """
        records = list(self._fetch(chromosome, start=start, end=end))
        return self._process_single_chromosome(chromosome, records)

    def fetch_regions(
        self, chromosome: str, regions: Iterable[Tuple[int, Optional[int]]]
    ) -> VariantTable:
        """
        Fetch records from a single chromosome that overlap the given regions.

        :param regions: a list of start, end tuples (end can be None)
        """
        records = []
        for start, end in regions:
            records.extend(list(self._fetch(chromosome, start=start, end=end)))
        return self._process_single_chromosome(chromosome, records)

    def __iter__(self) -> Iterator[VariantTable]:
        """
        Yield VariantTable objects for each chromosome.

        Multi-ALT sites are skipped. (TODO Have to fix that.)
        """
        ## self._vcf_reader is VariantFile object
        ## So it records is a list of VariantRecord objects
        for chromosome, records in itertools.groupby(self._vcf_reader, lambda record: record.chrom):
            yield self._process_single_chromosome(chromosome, records)

    def _process_single_chromosome(self, chromosome: str, records) -> VariantTable:
        n_snvs = 0
        n_other = 0
        n_multi = 0
        n_skip = 0  #To count the number of records that need to be skipped since they have more alleles than can be handled by Giggles
        table = VariantTable(chromosome, self.samples)
        prev_position = None
        ## records is a list of VariantRecord objects
        logger.info("Processing variants from Chromosome %s."%(chromosome))
        for record in records:
            if not record.alts:
                continue
            if len(record.alts) > 1:
                n_multi += 1
                
            pos, ref = record.start, str(record.ref)
            alts = record.alts
            if len(alts) > 15:
                n_skip += 1
                continue
            allele_origin = []
            for _, call in record.samples.items():
                allele_origin.append(call["GT"])
            allele_traversal = record.info["AT"]
            for alt in alts:
                if len(ref) == len(alt) == 1:
                    n_snvs += 1
                else:
                    n_other += 1

            if (prev_position is not None) and (prev_position > pos):
                raise VcfNotSortedError(
                    "VCF not ordered: {}:{} appears before {}:{}".format(
                        chromosome, prev_position + 1, chromosome, pos + 1
                    )
                )

            if prev_position == pos:
                warn_once(
                    logger, "Skipping duplicated position %s on chromosome %r", pos + 1, chromosome
                )
                continue
            prev_position = pos

            """
            Not reading phase, genotype or genotype quality information from the input vcf.
            Since this is not re-genotyping, the samples in the input vcf are not the ones
            we care about. We need info about the samples in the BAM/GAF file.
            Hence we are storing only None values for the number of samples (which is 1).
            """
            genotype_likelihoods = [None] * len(self.samples)
            genotypes = [Genotype([]) for i in range(len(self.samples))]
            phases = [None] * len(self.samples)
            
            variant = VcfVariant(position=pos, reference_allele=ref, alternative_allele=alts, allele_origin=allele_origin, allele_traversal=allele_traversal)
            table.add_variant(variant, genotypes, phases, genotype_likelihoods)

        logger.info("Processed Chromosome %s. Parsed %s SNVs, %s non-SNVs and %s multi-ALTs. Also skipped %s records exceeding max allele caparacity.", chromosome, n_snvs, n_other, n_multi, n_skip)

        return table

@dataclass
class VcfHeader:
    format_or_info: str
    id: str
    number: Union[str, int]
    typ: str
    description: str

    def line(self):
        return (
            "##{format_or_info}=<ID={id},Number={number},Type={typ},"
            'Description="{description}">'.format(
                format_or_info=self.format_or_info,
                id=self.id,
                number=self.number,
                typ=self.typ,
                description=self.description,
            )
        )


PREDEFINED_FORMATS = {
    "GL": VcfHeader(
        "FORMAT",
        "GL",
        "G",
        "Float",
        "Genotype Likelihood, log10-scaled likelihoods of the data given the"
        " called genotype for each possible genotype generated from the"
        " reference and alternate alleles given the sample ploidy",
    ),
    "GQ": VcfHeader("FORMAT", "GQ", 1, "Integer", "Phred-scaled genotype quality"),
    "GT": VcfHeader("FORMAT", "GT", 1, "String", "Genotype"),
    "HP": VcfHeader("FORMAT", "HP", ".", "String", "Phasing haplotype identifier"),
    "PQ": VcfHeader("FORMAT", "PQ", 1, "Float", "Phasing quality"),
    "PS": VcfHeader("FORMAT", "PS", 1, "Integer", "Phase set identifier"),
    "HS": VcfHeader("FORMAT", "HS", ".", "Integer", "Haploid phase set identifier"),
}

PREDEFINED_INFOS = {
    "AC": VcfHeader(
        "INFO",
        "AC",
        "A",
        "Integer",
        "Allele count in genotypes, for each ALT allele, in the same order as listed",
    ),
    "AN": VcfHeader("INFO", "AN", "A", "Integer", "Total number of alleles in called genotypes"),
    "END": VcfHeader("INFO", "END", 1, "Integer", "Stop position of the interval"),
    "SVLEN": VcfHeader(
        "INFO", "SVLEN", ".", "Integer", "Difference in length between REF and ALT alleles"
    ),
    "SVTYPE": VcfHeader("INFO", "SVTYPE", 1, "String", "Type of structural variant"),
}


def augment_header(header: VariantHeader, contigs: List[str], formats: List[str], infos: List[str]):
    """
    Add contigs, formats and infos to a VariantHeader.

    formats and infos are given as a list of strings, where each item is the ID of the header
    line to add. The full header info (Number, Type, Description) is taken from the PREDEFINED_*
    constants above. Any other FORMATs or INFOs that are not predefined will raise a VcfError.

    The header is modified in place.
    """
    for contig in contigs:
        header.contigs.add(contig)

    for fmt in formats:
        if fmt in header.formats:
            header.formats[fmt].remove_header()
        try:
            h = PREDEFINED_FORMATS[fmt]
        except KeyError:
            raise VcfError("FORMAT {!r} not defined in VCF header".format(fmt)) from None
        header.add_line(h.line())

    for info in infos:
        try:
            h = PREDEFINED_INFOS[info]
        except KeyError:
            raise VcfError("INFO {!r} not defined in VCF header".format(info)) from None
        header.add_line(h.line())


def missing_headers(path: str) -> Tuple[List[str], List[str], List[str]]:
    """
    Find contigs, FORMATs and INFOs that are used within the body of a VCF file, but are
    not listed in the header or that have an incorrect type.

    Return a tuple (contigs, formats, infos) where each of the items are lists of
    strings.

    The reason this function exists is that pysam.VariantFile crashes when we
    try to write a VCF record to it that uses contigs, INFOs or FORMATs that
    are missing from the header. See also
    <https://github.com/pysam-developers/pysam/issues/771>
    """
    with VariantFile(path) as variant_file:
        header = variant_file.header.copy()
        # Check for FORMATs that do not have the expected type
        incorrect_formats = []
        for fmt, v in variant_file.header.formats.items():
            if fmt not in PREDEFINED_FORMATS:
                continue
            h = PREDEFINED_FORMATS[fmt]
            if v.number != h.number or v.type != h.typ:
                if fmt == "PS" and v.type != h.typ:
                    raise VcfError(
                        "The input VCF/BCF contains phase set ('PS') tags that are of the"
                        " non-standard type '{}' instead of 'Integer'. Giggles cannot"
                        " overwrite these as it could produce inconsistent files."
                        " To proceed, you can use 'giggles unphase' to remove phasing"
                        " information from the input file".format(v.type)
                    )
                incorrect_formats.append(fmt)

        # Iterate through entire file and check which contigs, formats and
        # info fields are used
        contigs = []  # contigs encountered, in the proper order
        seen_contigs = set()
        formats = []  # FORMATs encountered, in the proper order
        seen_formats = set()
        seen_infos = set()  # INFOs encountered

        for record in variant_file:
            seen_infos.update(record.info)
            if record.alts is not None:
                for alt in record.alts:
                    # If there are "vague" ALT alleles such as <INS>, <DEL> etc, then
                    # the header needs to contain a LEN info entry even if LEN
                    # is never used
                    if alt.startswith("<"):
                        seen_infos.add("END")

            # For the contigs, we maintain a set *and* a list because we want to
            # keep track of the order of the contigs.
            if record.contig not in seen_contigs:
                contigs.append(record.contig)
            seen_contigs.add(record.contig)

            for fmt in record.format:
                if fmt not in seen_formats:
                    formats.append(fmt)
                seen_formats.add(fmt)

    # Determine which contigs are missing from the header
    header_contigs = set(header.contigs)
    missing_contigs = []
    for contig in contigs:
        if contig not in header_contigs:
            missing_contigs.append(contig)

    # Determine which FORMATs are missing from the header
    header_formats = set(header.formats)
    missing_formats = []
    for fmt in formats:
        if fmt in header_formats:
            continue
        missing_formats.append(fmt)

    # Determine which INFOs are missing from the header
    missing_infos = list(set(seen_infos) - set(header.info))

    return (missing_contigs, incorrect_formats + missing_formats, missing_infos)


@dataclass
class GenotypeChange:
    sample: str
    chromosome: str
    variant: VcfVariant
    old_gt: Genotype
    new_gt: Genotype


class VcfAugmenter(ABC):
    def __init__(
        self,
        in_path: str,
        bam_samples: Iterable[str],
        command_line: Optional[str],
        out_file: TextIO = sys.stdout,
        include_haploid_phase_sets: bool = False,
    ):
        """
        in_path -- Path to input VCF, used as template.
        command_line -- A string that will be added as a VCF header entry
            (use None to not add this to the VCF header)
        out_file -- Open file-like object to which VCF is written.
        tag -- which type of tag to write, either 'PS' or 'HP'. 'PS' is standardized;
            'HP' is compatible with GATKs ReadBackedPhasing.
        """
        # TODO This is slow because it reads in the entire VCF one extra time
        contigs, formats, infos = missing_headers(in_path)
        # TODO It would actually look nicer if the custom HS header was directly below PS
        if include_haploid_phase_sets and "HS" not in formats:
            formats.append("HS")
        # We repair the header (adding missing contigs, formats, infos) of the *input* VCF because
        # we will modify the records that we read, and these are associated with the input file.
        self._reader = VariantFile(in_path)
        augment_header(self._reader.header, contigs, formats, infos)
        if command_line is not None:
            command_line = '"' + command_line.replace('"', "") + '"'
            self._reader.header.add_meta("commandline", command_line)
        self._writer = VariantFile(out_file, mode="w", header=VariantHeader())
        self.setup_header(self._writer.header)
        for sample in bam_samples:
            self._writer.header.add_sample(sample)
        self._unprocessed_record: Optional[VariantRecord] = None
        self._reader_iter = iter(self._reader)

    @abstractmethod
    def setup_header(self, header):
        pass

    def close(self):
        self._writer.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    @property
    def samples(self) -> List[str]:
        return list(self._reader.header.samples)

    def _record_modifier(self, chromosome: str):
        for record in self._iterrecords(chromosome):
            new_record = self._writer.new_record(contig = record.contig, start = record.start, alleles = record.alleles, id = record.id, qual = record.qual, filter = record.filter, info = record.info)
            yield new_record
            self._writer.write(new_record)

    def _iterrecords(self, chromosome: str) -> Iterable[VariantRecord]:
        """Yield all records for the target chromosome"""
        n = 0
        if self._unprocessed_record is not None:
            assert self._unprocessed_record.chrom == chromosome
            yield self._unprocessed_record
            n += 1
        for record in self._reader_iter:
            n += 1
            if record.chrom != chromosome:
                # save it for later
                self._unprocessed_record = record
                assert n != 1
                return
            yield record


def genotype_code(gt: Optional[Tuple[Optional[int], ...]]) -> Genotype:
    """Return genotype encoded as PyVCF-compatible number"""
    if gt is None:
        result = Genotype([])
    elif any(allele is None for allele in gt):
        result = Genotype([])
    else:
        result = Genotype([allele for allele in gt])  # type: ignore
    return result


# class to print computed genotypes,likelihoods (still needs to be improved...)
# in input vcf, currently GT is still required..


class GenotypeVcfWriter(VcfAugmenter):
    """
    Read in a VCF file and write it back out with added genotyping information.

    Avoid reading in full chromosomes as that uses too much memory for
    multi-sample VCFs.
    """

    def __init__(self, in_path: str, bam_samples: Iterable[str], command_line: Optional[str], out_file: TextIO = sys.stdout):
        """
        in_path -- Path to input VCF, used as template.
        command_line -- A string that will be added as a VCF header entry.
        out_file -- Open file-like object to which VCF is written.
        """
        super().__init__(in_path, bam_samples, command_line, out_file)

    def setup_header(self, header: VariantHeader):
        """Called by baseclass constructor"""
        header.add_line(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype computed by Giggles genotyping algorithm">'
        )
        header.add_line(
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Phred-scaled genotype quality computed by Giggles genotyping algorithm">'
        )
        header.add_line(
            '##FORMAT=<ID=GL,Number=G,Type=Float,Description="Log10-scaled likelihoods for genotypes: 0/0, 0/1, 1/1, computed by Giggles genotyping algorithm">'
        )

        for fmt, v in self._reader.header.info.items():
            header.info.add(fmt,v.number, v.type, v.description)
        for contig in list(self._reader.header.contigs):
            header.contigs.add(contig)

    def write_genotypes(
        self, chromosome: str, variant_table: VariantTable, indels, ploidy: int = 2
    ) -> None:
        """
        Add genotyping information to all variants on a single chromosome.

        chromosome -- name of chromosome
        variant_table -- contains genotyping information for all accessible variant positions
        leave_unchanged -- if True, leaves records of current chromosome unchanged
        """

        # map positions to index
        genotyped_variants = dict()
        for i in range(len(variant_table)):
            genotyped_variants[variant_table.variants[i].position] = i

        # INT_TO_UNPHASED_GT = {0: (0, 0), 1: (0, 1), 2: (1, 1), -1: None}
        GT_GL_GQ = frozenset(["GT", "GL", "GQ"])
        for record in self._record_modifier(chromosome):
            pos = record.start
            if not record.alts:
                continue
            for sample, call in record.samples.items():
                geno = Genotype([])
                n_alleles = 1 + len(record.alts)
                n_genotypes = binomial_coefficient(ploidy + n_alleles - 1, n_alleles - 1)
                geno_l = [1 / n_genotypes] * int(n_genotypes)
                geno_q = None

                # for genotyped variants, get computed likelihoods/genotypes (for all others, give uniform likelihoods)
                if pos in genotyped_variants:
                    likelihoods = variant_table.genotype_likelihoods_of(sample)[
                        genotyped_variants[pos]
                    ]
                    # likelihoods can be 'None' if position was not accessible
                    if likelihoods is not None:
                        geno_l = [l for l in likelihoods]  # type: ignore
                        geno = variant_table.genotypes_of(sample)[genotyped_variants[pos]]

                # Compute GQ
                geno_index = geno.get_index()
                geno_q = sum(geno_l[i] for i in range(n_genotypes) if i != geno_index)
                # TODO default value ok?
                # store likelihoods log10-scaled

                # Temporarily overwrite the GT field with a (fake) genotype that indicates a
                # diploid sample. Otherwise, if the GT field happens to be empty, pysam
                # complains that we are setting an incorrect number of GL values.
                call["GT"] = tuple([0] * ploidy)

                call["GL"] = [max(math.log10(j), -1000) if j > 0 else -1000 for j in geno_l]
                call["GT"] = tuple(geno.as_vector())

                # store quality as phred score
                if not geno.is_none():
                    # TODO default value ok?
                    assert geno_q is not None
                    if geno_q > 0:
                        call["GQ"] = min(round(-10.0 * math.log10(geno_q)), 10000)
                    else:
                        call["GQ"] = 10000
                else:
                    call["GQ"] = None

                record.qual = None

                # delete all other genotype information that might have been present before
                for tag in set(call.keys()) - GT_GL_GQ:
                    del call[tag]
