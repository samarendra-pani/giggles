"""
Functions for reading VCFs.
"""
# Code modified from WhatsHap (https://github.com/whatshap/whatshap)

import os
import sys
import math
import itertools
from dataclasses import dataclass
from abc import ABC, abstractmethod
from os import PathLike
from typing import List, Sequence, Tuple, Iterable, Optional, Union, TextIO, Iterator

from pysam import VariantFile, VariantHeader, VariantRecord

from .core import (
    GenotypeLikelihoods,
    Genotype,
    binomial_coefficient,
    get_max_genotype_alleles
)
from .align import edit_distance

from .logger import logger, warn_once

@dataclass
class VariantCallPhase:
    block_id: int  # numeric id of the phased block
    phase: Tuple[int, ...]  # alleles representing the phasing. (1, 0) is 1|0
    quality: Optional[int]


class VcfVariant:
    """A variant in a VCF file (not to be confused with core.Variant)"""

    __slots__ = ("id", "position", "position_on_ref", "reference_allele", "alternative_allele", "allele_origin", "length_on_path", "state", "distance_matrix")

    def __init__(self, id: str, position: int, reference_allele: str, alternative_allele: tuple, allele_origin: list, use_distance_matrix: bool = False):
        
        self.id = id
        # This is the position on the backbone reference (the position given in the VCF in the 0-base)
        # The position will be 1 less than what is seen in the VCF
        self.position_on_ref = position
        # This is the position on the paths (the position used to find the variant locations on paths).
        # This changes for every new alignment path.
        self.position = None
        self.reference_allele = reference_allele    # reference allele given in the VCF
        self.alternative_allele = alternative_allele    # alternate alleles given in the VCF
        self.allele_origin = allele_origin  # the phased genotypes in the VCF. Used for the HMM
        # This is the length of the variant on the alignment path.
        # This is needed since the CIGAR string processing needs this length.
        # This changes for every new alignment.
        self.length_on_path = None
        # the 'state' variable is used to denote how this variant is covered by the read.
        # this changes for every new alignment
        # the following values are possible:
        #   - 0: the whole variant is covered.
        #   - 1: the read starts within this variant.
        #   - 2: the read ends within this variant.
        #   - 3: the read starts and ends within this variant.
        self.state = None
        if use_distance_matrix and self.is_sv():
            self.distance_matrix = self.calculate_distance_matrix()
        else:
            self.distance_matrix = None

    #def __repr__(self):
    #    return "VcfVariant({}, {}, {}, {!r}, {!r}, {!r})".format(
    #        self.id, self.position_on_ref, self.position_on_ref + len(self.reference_allele), self.reference_allele, self.alternative_allele, self.allele_origin
    #    )
    
    def __repr__(self):
        return "VcfVariant({}, {}, {})".format(
            self.id, self.position_on_ref, self.position_on_ref + len(self.reference_allele)
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
    
    def is_sv(self) -> bool:
        return True if any(len(alt) > 50 or len(self.reference_allele) > 50 for alt in self.alternative_allele) else False

    def get_variant_bo(self, rgfa):
        if not self.is_sv():
            # no tag available for external variants
            return None
        start_id = self.id.split('>')[1]
        end_id = self.id.split('>')[-1]
        start_bo = rgfa.get_node(start_id).tags['BO']
        end_bo = rgfa.get_node(end_id).tags['BO']
        assert end_bo == start_bo + 2, f"Inconsistent BO tags for bubble {self.id}"
        return start_bo + 1

    def has_anchor_base(self):
        anchor_base = self.reference_allele[0]
        has_anchor = True
        for alt in self.alternative_allele:
            if alt[0] != anchor_base:
                has_anchor = False
                break
        return has_anchor
    
    def remove_anchor_base(self):
        # removing the anchor base that was added in the VCF
        self.position_on_ref += 1
        self.reference_allele = self.reference_allele[1:]
        new_alts = []
        for alt in self.alternative_allele:
            new_alts.append(alt[1:])
        self.alternative_allele = tuple(new_alts)
    
    # Calculate distance estimates between alleles. Max distance of 30 is considered.
    # storing the distance in a 1D array using canonical index
    def calculate_distance_matrix(self):
        self.distance_matrix = []
        alleles = [self.reference_allele] + self.alternative_allele
        n = len(alleles)
        for i in range(len(alleles)):
            for j in range(i+1, len(alleles)):
                k = (i * ((2*n) - i - 1))/2 + j - i - 1
                assert (len(self.distance_matrix) == k)     # checking for correctness of cannonical index.
                allele1 = alleles[i]
                allele2 = alleles[j]
                if abs(len(allele1) - len(allele2)) >= 50:
                    self.distance_matrix.push(50)
                    continue
                self.distance_matrix.push(edit_distance(allele1, allele2, 50))
    
    # get the distance between allele i and allele j
    # converting i and j into canonical index
    def get_distance(self, i, j):
        if i == j:
            return 0
        assert i < j
        n = len(self.alternative_allele) + 1
        k = (i * ((2*n) - i - 1))/2 + j - i - 1
        return self.distance_matrix[k]


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
        self.variants = []
        
        # Separate lists for VCF samples and GAF/BAM sample
        self._sample_to_index = {sample: index for index, sample in enumerate(samples)}

        self.query_genotypes= []
        self.query_genotype_likelihoods = []

    def __len__(self) -> int:
        return len(self.variants)

    def add_variant(
        self,
        variant: VcfVariant,
    ) -> None:
        self.variants.append(variant)
        
        # Adding empty Genotype object for the GAF/BAM sample
        self.query_genotypes.append(Genotype([]))
        self.query_genotype_likelihoods.append(None)

    def id_of(self, sample: str) -> int:
        """Return a unique int id of a sample given by name"""
        return self._sample_to_index[sample]
    

    # Making a copy of all the functions to include GAF/BAM sample functions separately.
    def query_genotypes_of(self) -> List[Genotype]:
        """Retrieve genotypes by sample name"""
        return self.query_genotypes

    def query_set_genotypes_of(self, genotypes: List[Genotype]) -> None:
        """Set genotypes by sample name"""
        assert len(genotypes) == len(self.variants)
        self.query_genotypes = genotypes

    def query_genotype_likelihoods_of(self) -> List[Optional[GenotypeLikelihoods]]:
        """Retrieve genotype likelihoods by sample name"""
        return self.query_genotype_likelihoods

    def query_set_genotype_likelihoods_of(
        self, genotype_likelihoods: List[Optional[GenotypeLikelihoods]]
    ) -> None:
        """Set genotype likelihoods by sample name"""
        assert len(genotype_likelihoods) == len(self.variants)
        self.query_genotype_likelihoods = genotype_likelihoods


# TODO: Clean up VcfReader.
# There is space for haplotagging but there is no plans for developing in that direction yet.
class VcfReader:
    """
    Read a VCF file chromosome by chromosome.
    """

    def __init__(
        self,
        path: Union[str, PathLike],
        indels: bool = False,
        required_chr: List = None,
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
        self.vcf_samples = list(self._vcf_reader.header.samples)
        self.required_chr = required_chr
        
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
                raise Exception("Invalid chromosome found in VCF file")
            elif "fetch requires an index" in e.args[0]:
                raise Exception(f"{self._path} is missing an index (.tbi or .csi)")
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
        # self._vcf_reader is pysam.VariantFile object
        # So it records is a list of VariantRecord objects
        # Problem with itertools.groupby is in its processing of the entire VCF file instead of being a generator yielding records.
        # Requires high memory initially
        # TODO: Possible to make this more memory efficient?
        for chromosome, records in itertools.groupby(self._vcf_reader, lambda record: record.chrom):
            if (not self.required_chr) or (chromosome in self.required_chr):
                logger.info(f"======== Working on chromosome {chromosome}")
                yield self._process_single_chromosome(chromosome, records)
            else:
                logger.info(f"======== Skipping chromosome {chromosome}")

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
        n_snvs = 0
        n_other = 0
        n_multi = 0
        n_ext = 0
        n_skip = 0  #To count the number of records that need to be skipped since they have more alleles than can be handled by Giggles
        table = VariantTable(chromosome, self.vcf_samples)
        prev_position = None
        ## records is a list of VariantRecord objects
        logger.info(f"Processing variants from Chromosome {chromosome}.")
        for record in records:
            if not record.alts:
                continue
            id = record.id
            if len(record.alts) > 1:
                n_multi += 1
                
            pos, ref = record.start, str(record.ref)
            alts = record.alts
            if len(alts) >= get_max_genotype_alleles():
                # logger.warning(f'Skipping position {pos} of chromosome {chromosome}. Position has more alleles than currently supported.')
                n_skip += 1
                continue
            allele_origin = []
            for _, call in record.samples.items():
                allele_origin.append(call["GT"])
            if id.__contains__('EXT'):
                n_ext += 1
            for alt in alts:
                if len(ref) == len(alt) == 1:
                    n_snvs += 1
                else:
                    n_other += 1

            if (prev_position is not None) and (prev_position > pos):
                raise Exception(
                    f"VCF not ordered: {chromosome}:{prev_position + 1} appears before {chromosome}:{pos + 1}"
                )

            if prev_position == pos:
                warn_once(
                    logger, "Skipping duplicated position %s on chromosome %r", pos + 1, chromosome
                )
                continue
            prev_position = pos
            
            variant = VcfVariant(id = id, position=pos, reference_allele=ref, alternative_allele=alts, allele_origin=allele_origin)
            if variant.has_anchor_base():
                variant.remove_anchor_base()
            table.add_variant(variant)

        logger.info(f"Processed Chromosome {chromosome}. Parsed {n_snvs} SNVs, {n_other} non-SNVs and {n_multi} multi-ALTs. Identified {n_ext} external variants added to the graph variants. Also skipped {n_skip} records exceeding max allele caparacity.")

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
            raise Exception(f"FORMAT {fmt} not defined in VCF header")
        header.add_line(h.line())

    for info in infos:
        try:
            h = PREDEFINED_INFOS[info]
        except KeyError:
            raise Exception(f"INFO {info} not defined in VCF header")
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
                    raise Exception(
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
        sample: str,
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

    def __init__(self, in_path: str, sample: str, command_line: Optional[str], out_file: TextIO = sys.stdout):
        """
        in_path -- Path to input VCF, used as template.
        command_line -- A string that will be added as a VCF header entry.
        out_file -- Open file-like object to which VCF is written.
        """
        super().__init__(in_path, sample, command_line, out_file)

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
                    likelihoods = variant_table.query_genotype_likelihoods_of()[
                        genotyped_variants[pos]
                    ]
                    # likelihoods can be 'None' if position was not accessible
                    if likelihoods is not None:
                        geno_l = [l for l in likelihoods]  # type: ignore
                        geno = variant_table.query_genotypes_of()[genotyped_variants[pos]]

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
