"""
Functions for reading VCFs.
"""
# Code modified from WhatsHap (https://github.com/whatshap/whatshap)

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
    get_max_genotype_ploidy,
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

    __slots__ = ("id", "position", "position_on_ref", "reference_allele", "alternative_allele", "allele_origin", "allele_traversal","length_on_path")

    def __init__(self, id: str, position: int, reference_allele: str, alternative_allele: tuple, allele_origin: list, allele_traversal: tuple):
        
        self.id = id
        self.position_on_ref = position     # This is the position on the backbone reference (the position given in the VCF)
        self.position = position            # This is the position on the paths (the position used to find the variant locations on paths). This changes for every new alignment path.
        self.reference_allele = reference_allele
        self.alternative_allele = alternative_allele
        self.allele_origin = allele_origin
        self.allele_traversal = allele_traversal
        self.length_on_path = None
        
    def __repr__(self):
        return "VcfVariant({}, {}, {!r}, {!r}, {!r})".format(
            self.id, self.position_on_ref, self.reference_allele, self.alternative_allele, self.allele_origin
        )

    def __hash__(self):
        return hash((self.position_on_ref, self.reference_allele, self.alternative_allele))

    def __eq__(self, other):
        return (
            (self.position_on_ref == other.position_on_ref)
            and (self.reference_allele == other.reference_allele)
            and (self.alternative_allele == other.alternative_allele)
        )

    def __lt__(self, other):
        return (self.position_on_ref, self.reference_allele, self.alternative_allele) < (
            other.position_on_ref,
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

    def __init__(self, chromosome: str, query_samples: List[str], samples: List[str]):
        self.chromosome = chromosome
        self.samples = samples
        self.query_samples = query_samples
        self.variants: List[VcfVariant] = []
        
        # Separate lists for VCF samples and GAF/BAM sample
        self.genotypes: List[List[Genotype]] = [[] for _ in samples]
        self.phases: List[List[Optional[VariantCallPhase]]] = [[] for _ in samples]
        self.genotype_likelihoods: List[List[Optional[GenotypeLikelihoods]]] = [[] for _ in samples]
        self._sample_to_index = {sample: index for index, sample in enumerate(samples)}

        self._query_sample_to_index = {sample: index for index, sample in enumerate(query_samples)}
        self.query_genotypes: List[List[Genotype]] = [[] for _ in query_samples]
        self.query_genotype_likelihoods: List[List[Optional[GenotypeLikelihoods]]] = [[] for _ in query_samples]

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
        
        # Adding info for the vcf samples
        for i, genotype in enumerate(genotypes):
            assert isinstance(genotype, Genotype)
            self.genotypes[i].append(genotype)
        for i, phase in enumerate(phases):
            self.phases[i].append(phase)
        for i, gl in enumerate(genotype_likelihoods):
            self.genotype_likelihoods[i].append(gl)
        
        # Adding empty Genotype object for the GAF/BAM sample
        for i in range(len(self.query_samples)):
            self.query_genotypes[i].append(Genotype([]))
        for i in range(len(self.query_samples)):
            self.query_genotype_likelihoods[i].append(None)

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
    

    # Making a copy of all the functions to include GAF/BAM sample functions separately.
    def query_genotypes_of(self, sample: str) -> List[Genotype]:
        """Retrieve genotypes by sample name"""
        return self.query_genotypes[self._query_sample_to_index[sample]]

    def query_set_genotypes_of(self, sample: str, genotypes: List[Genotype]) -> None:
        """Set genotypes by sample name"""
        assert len(genotypes) == len(self.variants)
        self.query_genotypes[self._query_sample_to_index[sample]] = genotypes

    def query_genotype_likelihoods_of(self, sample: str) -> List[Optional[GenotypeLikelihoods]]:
        """Retrieve genotype likelihoods by sample name"""
        return self.query_genotype_likelihoods[self._query_sample_to_index[sample]]

    def query_set_genotype_likelihoods_of(
        self, sample: str, genotype_likelihoods: List[Optional[GenotypeLikelihoods]]
    ) -> None:
        """Set genotype likelihoods by sample name"""
        assert len(genotype_likelihoods) == len(self.variants)
        self.query_genotype_likelihoods[self._query_sample_to_index[sample]] = genotype_likelihoods

    def query_id_of(self, sample: str) -> int:
        """Return a unique int id of a sample given by name"""
        return self._query_sample_to_index[sample]


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
        ignore_genotypes: bool = True,
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
        self.samples = bam_samples 
        self.vcf_samples = list(self._vcf_reader.header.samples)
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

    @staticmethod
    def _extract_HP_phase(call) -> Optional[VariantCallPhase]:
        hp = call.get("HP")
        if hp is None or hp == (".",):
            return None
        fields = [[int(x) for x in s.split("-")] for s in hp]
        for i in range(len(fields)):
            assert fields[0][0] == fields[i][0]
        block_id = fields[0][0]
        order = [field[1] - 1 for field in fields]
        phase = call["GT"]
        phase = tuple(phase[order.index(i)] for i in range(len(order)))
        return VariantCallPhase(block_id=block_id, phase=phase, quality=call.get("PQ", None))

    @staticmethod
    def _extract_GT_PS_phase(call) -> Optional[VariantCallPhase]:
        if not call.phased:
            return None
        is_het = not all(x == call["GT"][0] for x in call["GT"])
        if not is_het:
            return None
        block_id = call.get("PS", 0)
        phase = call["GT"]
        return VariantCallPhase(block_id=block_id, phase=phase, quality=call.get("PQ", None))

    def _process_single_chromosome(self, chromosome: str, records) -> VariantTable:
        phase_detected = None
        n_snvs = 0
        n_other = 0
        n_multi = 0
        n_skip = 0  #To count the number of records that need to be skipped since they have more alleles than can be handled by Giggles
        table = VariantTable(chromosome, self.samples, self.vcf_samples)
        prev_position = None
        ## records is a list of VariantRecord objects
        logger.info("Processing variants from Chromosome %s."%(chromosome))
        for record in records:
            if not record.alts:
                continue
            id = record.id
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
            For Genotyping:
            Not reading phase, genotype or genotype quality information from the input vcf.
            Since this is not re-genotyping, the samples in the input vcf are not the ones
            we care about. We need info about the samples in the BAM/GAF file.
            Hence we are storing only None values for the number of samples (which is 1).

            For Haplotagging:
            Need the phase and genotype information since the haplotagging needs that.
            """
            # Read phasing information (allow GT/PS or HP phase information, but not both),
            # if requested
            if self._phases:
                phases = []
                for call in record.samples.values():
                    phase = None
                    for extract_phase, phase_name in [
                        (self._extract_HP_phase, "HP"),
                        (self._extract_GT_PS_phase, "GT_PS"),
                    ]:
                        p = extract_phase(call)
                        if p is not None:
                            if phase_detected is None:
                                phase_detected = phase_name
                            elif phase_detected != phase_name:
                                raise MixedPhasingError(
                                    "Mixed phasing information in input VCF (e.g. mixing PS "
                                    "and HP fields)"
                                )
                            phase = p
                            # check for ploidy consistency and limits
                            phase_ploidy = len(p.phase)
                            if phase_ploidy > get_max_genotype_ploidy():
                                raise PloidyError(
                                    "Ploidies higher than {} are not supported."
                                    "".format(get_max_genotype_ploidy())
                                )
                            elif p is None or p.block_id is None or p.phase is None:
                                pass
                            elif self.ploidy is None:
                                self.ploidy = phase_ploidy
                            elif phase_ploidy != self.ploidy:
                                print(f"phase= {phase}")
                                raise PloidyError(
                                    "Phasing information contains inconsistent ploidy ({} and "
                                    "{})".format(self.ploidy, phase_ploidy)
                                )
                    phases.append(phase)
            else:
                phases = [None] * len(record.samples)
            
            if not self._ignore_genotypes:
                # check for ploidy consistency and limits
                genotype_lists = [call.get("GT", None) for call in record.samples.values()]
                for geno in genotype_lists:
                    if geno is None or None in geno:
                        continue
                    geno_ploidy = len(geno)
                    if geno_ploidy > get_max_genotype_ploidy():
                        raise PloidyError(
                            "Ploidies higher than {} are not supported."
                            "".format(get_max_genotype_ploidy())
                        )
                    elif self.ploidy is None:
                        self.ploidy = geno_ploidy
                    elif geno_ploidy != self.ploidy:
                        raise PloidyError(
                            "Inconsistent ploidy ({} and " "{})".format(self.ploidy, geno_ploidy)
                        )

                genotypes = [genotype_code(geno_list) for geno_list in genotype_lists]
            else:
                genotypes = [Genotype([]) for _ in self.vcf_samples]
                phases = [None] * len(self.vcf_samples)
            
            genotype_likelihoods = [None] * len(self.vcf_samples)
            
            variant = VcfVariant(id = id, position=pos, reference_allele=ref, alternative_allele=alts, allele_origin=allele_origin, allele_traversal=allele_traversal)
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
            new_record = self._writer.new_record(contig = record.contig, 
                                                 start = record.start, 
                                                 alleles = record.alleles, 
                                                 id = record.id, 
                                                 qual = record.qual, 
                                                 filter = record.filter, 
                                                 info = record.info)
            yield new_record
            self._writer.write(new_record)

    def _iterrecords(self, chromosome: str) -> Iterable[VariantRecord]:
        """Yield all records for the target chromosome"""
        n = 0
        for record in self._reader.fetch(contig=chromosome):
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
            genotyped_variants[variant_table.variants[i].position_on_ref] = i

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
                    likelihoods = variant_table.query_genotype_likelihoods_of(sample)[
                        genotyped_variants[pos]
                    ]
                    # likelihoods can be 'None' if position was not accessible
                    if likelihoods is not None:
                        geno_l = [l for l in likelihoods]  # type: ignore
                        geno = variant_table.query_genotypes_of(sample)[genotyped_variants[pos]]

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
