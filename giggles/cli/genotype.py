"""
Genotype variants

Runs only the genotyping algorithm. Genotype Likelihoods are computed using the
forward backward algorithm.
"""
import logging
import sys
import platform
from typing import Sequence
from collections import defaultdict

from contextlib import ExitStack

from giggles import __version__
from giggles import vcf
from giggles.vcf import VcfReader, GenotypeVcfWriter
from giggles.core import (
    GenotypeHMM,
    ReadSet,
    Pedigree,
    NumericSampleIds,
    PhredGenotypeLikelihoods,
    Genotype
)
from giggles.pedigree import (
    UniformRecombinationCostComputer,
)
from giggles.timer import StageTimer
from giggles.cli import log_memory_usage
from giggles.phase import select_reads, setup_families
from giggles.cli import CommandLineError, PhasedInputReader, read_haplotags


logger = logging.getLogger(__name__)

def bin_coeff(n, k):
    if (k < 0) or (n < 0) or (n < k):
        return 0
	
    result = 1.0
    if (k > n-k):
        k = n-k
	
    for i in range(k):
        result *= (n-i)
        result /= (i+1)
	
    return int(result)


def int_to_diploid_multiallelic_gt(numeric_repr):
    """Converts the classic numeric representation of multi-allelic, diploid genotypes
    into a genotype object
    """
    if numeric_repr == -1:
        return Genotype([])
    ploidy = 2
    genotype = [-1,-1]
    pth = ploidy
    max_allele_index = numeric_repr
    leftover_genotype_index = numeric_repr

    while (pth > 0):
        for allele_index in range(max_allele_index+1):
            i = bin_coeff(pth + allele_index - 1, pth)
            if (i >= leftover_genotype_index) or (allele_index == max_allele_index):
                if (i > leftover_genotype_index):
                    allele_index -= 1
                leftover_genotype_index -= bin_coeff(pth + allele_index - 1, pth)
                pth -= 1
                max_allele_index = allele_index
                genotype[pth] = allele_index
                break
    
    return Genotype(genotype)


def determine_genotype(likelihoods: Sequence[float], threshold_prob: float, n_allele: int) -> float:
    """given genotype likelihoods for 0/0, 0/1, 1/1, determines likeliest genotype"""

    assert bin_coeff(n_allele + 1, n_allele - 1) == len(likelihoods)
    to_sort = []
    for i in range(len(likelihoods)):
        to_sort.append((likelihoods[int_to_diploid_multiallelic_gt(i)], i))
    to_sort.sort(key=lambda x: x[0])

    # make sure there is a unique maximum which is greater than the threshold
    if (to_sort[-1][0] > to_sort[-2][0]) and (to_sort[-1][0] > threshold_prob):
        return int_to_diploid_multiallelic_gt(to_sort[-1][1])
    else:
        return int_to_diploid_multiallelic_gt(-1)


def run_genotype(
    mapped_read_files,
    variant_file,
    reference_fasta=None,
    rgfa=None,
    read_fasta=None,
    haplotag_tsv=None,
    keep_untagged=False,
    output=sys.stdout,
    samples=None,
    chromosomes=None,
    mapping_quality=20,
    max_coverage=15,
    recombrate=1.26,
    gt_qual_threshold=0,
    realign_mode="ed",
    overhang=10,
    gap_start=3,
    gap_extend=1,
    mismatch=2,
    match_probability=0.85,
    mismatch_probability=0.05,
    insertion_probability=0.05,
    deletion_probability=0.05,
    write_command_line_header=True,
    eff_pop_size = 10
):
    """
    For now: this function only runs the genotyping algorithm. Genotype likelihoods for
    all variants are computed using the forward backward algorithm
    """
    timers = StageTimer()
    logger.info(
        "This is Giggles (genotyping) %s running under Python %s",
        __version__,
        platform.python_version(),
    )
    if write_command_line_header:
        command_line = "(giggles {}) {}".format(__version__, " ".join(sys.argv[1:]))
    else:
        command_line = None
    with ExitStack() as stack:
        # read the given input files (BAMs, VCFs, ref...)
        numeric_sample_ids = NumericSampleIds()
        phased_input_reader = stack.enter_context(
            PhasedInputReader(
                mapped_read_files,
                reference_fasta,
                rgfa,
                read_fasta,
                numeric_sample_ids,
                mapq_threshold=mapping_quality,
                realign_mode=realign_mode,
                overhang=overhang,
                gap_start=gap_start,
                gap_extend=gap_extend,
                default_mismatch=mismatch,
                em_prob_params=[match_probability, mismatch_probability, insertion_probability, deletion_probability]
            )
        )

        haplotags = read_haplotags(haplotag_tsv)
        logger.debug("Initial Parsing of Alignments Done. PhasedInputReader object successfully created.")

        # vcf writer for final genotype likelihoods
        vcf_writer = stack.enter_context(GenotypeVcfWriter(command_line=command_line, in_path=variant_file, out_file=output, bam_samples = samples))
        
        # The samples in the gaf or bam is given as input since it will be used to make variant tables with those samples.
        # The variant tables are then simply just updated after the HMM is run.
        vcf_reader = stack.enter_context(
            VcfReader(
                variant_file, bam_samples=samples, indels=True, genotype_likelihoods=False, phases=False
            )
        )

        recombination_cost_computer = UniformRecombinationCostComputer(recombrate, eff_pop_size)

        samples = frozenset(samples)
        families, family_trios = setup_families(samples, max_coverage)
        for trios in family_trios.values():
            for trio in trios:
                # Ensure that all mentioned individuals have a numeric id
                _ = numeric_sample_ids[trio.child]

        # compute genotype likelihood threshold
        gt_prob = 1.0 - (10 ** (-gt_qual_threshold / 10.0))

        # Count number of samples are present in the multisample reference graph vcf file 
        n_samples = len(list(vcf_reader._vcf_reader.header.samples))
        
        # Iterating over chromosomes present in the multisample reference graph vcf file
        for variant_table in timers.iterate("parse_vcf", vcf_reader):

            chromosome = variant_table.chromosome
            if (not chromosomes) or (chromosome in chromosomes):
                logger.info("======== Working on chromosome %r", chromosome)
            else:
                continue
            
            # create a mapping of genome positions to indices
            var_pos_to_ind = dict()
            n_allele_position = dict()
            allele_references = dict()
            ids = dict()
            for i in range(len(variant_table.variants)):
                var_pos_to_ind[variant_table.variants[i].position_on_ref] = i
                v = variant_table.variants[i]
                n_allele_position[v.position_on_ref] = len(v.alternative_allele)+1      ##Contains the number of alleles at every variant position
                allele_references[v.position_on_ref] = v.allele_origin
                ids[v.position_on_ref] = v.id            
            #Prior genotyping with equal probabilities
            for sample in samples:
                variant_table.query_set_genotype_likelihoods_of(
                    sample, [PhredGenotypeLikelihoods([1/(bin_coeff(n_allele_position[pos] + 1, n_allele_position[pos] - 1))] * (bin_coeff(n_allele_position[pos] + 1, n_allele_position[pos] - 1)) , 2, n_allele_position[pos]) for pos in list(var_pos_to_ind.keys())]
                )
            
            # Iterate over all families to process, i.e. a separate DP table is created
            # for each family.
            for representative_sample, family in sorted(families.items()):
                if len(family) == 1:
                    logger.info("---- Processing individual %s", representative_sample)
                else:
                    logger.info("---- Processing family with individuals: %s", ",".join(family))
                max_coverage_per_sample = max(1, max_coverage // len(family))
                logger.info("Using maximum coverage per sample of %dX", max_coverage_per_sample)
                trios = family_trios[representative_sample]
                assert (len(family) == 1) or (len(trios) > 0)

                # Get the reads belonging to each sample
                readsets = dict()
                for sample in family:
                    with timers("read_alignment"):
                        readset = phased_input_reader.read(
                            chromosome, variant_table.variants, sample, haplotags, keep_untagged
                        )

                    with timers("select"):
                        readset = readset.subset(
                            [i for i, read in enumerate(readset) if len(read) >= 2]
                        )
                        logger.info(
                            "Kept %d reads that cover at least two variants each", len(readset)
                        )
                        selected_reads = select_reads(readset, max_coverage_per_sample)
                    readsets[sample] = selected_reads
                    
                # Merge reads into one ReadSet (note that each Read object
                # knows the sample it originated from).
                all_reads = ReadSet()
                for sample, readset in readsets.items():
                    for read in readset:
                        if not read.is_sorted():
                            read.sort()
                        all_reads.add(read)

                all_reads.sort()
                
                # Determine which variants can (in principle) be phased
                accessible_positions = sorted(all_reads.get_positions())
                accessible_positions_n_allele = []
                accessible_positions_allele_references = []
                for index, position in enumerate(accessible_positions):
                    accessible_positions_n_allele.append(n_allele_position[position])
                    allele_reference_to_list = []
                    for ref_sample in allele_references[position]:
                        for hap in ref_sample:
                            try:
                                allele_reference_to_list.append(int(hap))
                            except TypeError:
                                allele_reference_to_list.append(-1)
                    accessible_positions_allele_references.append(allele_reference_to_list)
                logger.info(
                    "Variants covered by at least one phase-informative "
                    "read in at least one individual after read selection: %d",
                    len(accessible_positions),
                )
                
                # Create Pedigree
                pedigree = Pedigree(numeric_sample_ids)
                for sample in family:
                    # genotypes are assumed to be unknown, so ignore information that
                    # might already be present in the input vcf
                    all_genotype_likelihoods = variant_table.query_genotype_likelihoods_of(sample)
                    genotype_l = [
                        all_genotype_likelihoods[var_pos_to_ind[a_p]] for a_p in accessible_positions
                    ]
                    pedigree.add_individual(
                        sample, [Genotype([]) for i in range(len(accessible_positions))], genotype_l
                    )
                for trio in trios:
                    pedigree.add_relationship(
                        father_id=trio.father, mother_id=trio.mother, child_id=trio.child
                    )

                recombination_costs = recombination_cost_computer.compute(accessible_positions)
                # Finally, run genotyping algorithm
                with timers("genotyping"):
                    problem_name = "genotyping"
                    logger.info(
                        "Genotype %d sample%s by solving the %s problem ...",
                        len(family),
                        "s" if len(family) > 1 else "",
                        problem_name,
                    )
                    forward_backward_table = GenotypeHMM(
                        numeric_sample_ids,
                        all_reads,
                        recombination_costs,
                        pedigree,
                        2*n_samples,
                        accessible_positions,
                        accessible_positions_n_allele,
                        accessible_positions_allele_references
                    )
                    
                    # store results
                    for s in family:
                        likelihood_list = variant_table.query_genotype_likelihoods_of(s)
                        genotypes_list = variant_table.query_genotypes_of(s)

                        for pos in range(len(accessible_positions)):
                            likelihoods = forward_backward_table.get_genotype_likelihoods(s, pos, accessible_positions_n_allele[pos])
                            # compute genotypes from likelihoods and store information
                            geno = determine_genotype(likelihoods, gt_prob, accessible_positions_n_allele[pos])
                            assert isinstance(geno, Genotype)
                            genotypes_list[var_pos_to_ind[accessible_positions[pos]]] = geno
                            likelihood_list[var_pos_to_ind[accessible_positions[pos]]] = likelihoods

                        variant_table.query_set_genotypes_of(s, genotypes_list)
                        variant_table.query_set_genotype_likelihoods_of(s, likelihood_list)

            with timers("write_vcf"):
                logger.info("======== Writing VCF")
                vcf_writer.write_genotypes(chromosome, variant_table, indels=True)
                logger.info("Done writing VCF")

            logger.debug("Chromosome %r finished", chromosome)

    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info("Time spent reading alignments:               %9.2f s", timers.elapsed("read_alignment"))
    logger.info("Time spent parsing VCF:                      %9.2f s", timers.elapsed("parse_vcf"))
    logger.info("Time spent selecting reads:                  %9.2f s", timers.elapsed("select"))
    logger.info(
        "Time spent genotyping:                          %9.2f s", timers.elapsed("genotyping")
    )
    logger.info("Time spent writing VCF:                      %9.2f s", timers.elapsed("write_vcf"))
    logger.info("Time spent on rest:                          %9.2f s", total_time - timers.sum())
    logger.info("Total elapsed time:                          %9.2f s", total_time)
    logger.info("Total elapsed time:                          %9.2f hr", total_time/3600)


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('variant_file', metavar='VCF', help='VCF file with variants to be genotyped (can be gzip-compressed)')
    arg('mapped_read_files', nargs='*', metavar='READS', help='BAM or GAF file. For GAF file, it expects an index created by scaffold-sort.py with .gai extension.')

    arg('-o', '--output', default=sys.stdout,
        help='Output VCF file. Add .gz to the file name to get compressed output. '
        'If omitted, use standard output.')
    arg('--reference-fasta', '-rf', metavar='FASTA',
        help='FASTA Reference file. Provide this with the BAM file. '
        'If no index (.fai) exists, it will be created')
    arg('--rgfa', metavar='rGFA',
        help='GFA Reference file. Provide this along with the input GAF file')
    arg('--read-fasta', '-f', metavar='FASTA',
        help='FASTA file with the reads used in GAF file. Provide this with the GAF file unless the read sequences are given in the GAF file with RS or rs tag. '
        'If no index (.fai) exists, it will be created')
    arg('--haplotag-tsv', metavar='HAPLOTAG', help='TSV file containing the haplotag and phaseset information.')

    arg = parser.add_argument_group('Input pre-processing, selection and filtering').add_argument
    arg('--max-coverage', '-H', metavar='MAX_COV', default=15, type=int,
        help='Reduce coverage to at most MAX_COV (default: %(default)s).')
    arg('--mapping-quality', '--mapq', metavar='QUAL',
        default=20, type=int, help='Minimum mapping quality (default: %(default)s)')
    arg('--sample', dest='samples', metavar='SAMPLE', default=[], action='append',
        help='Name of a sample to genotype. Has to be given.')
    arg('--chromosome', dest='chromosomes', metavar='CHROMOSOME', default=[], action='append',
        help='Name of chromosome to genotyped. If not given, all chromosomes in the '
        'input VCF are genotyped. Can be used multiple times.')
    arg('--gt-qual-threshold', metavar='GT_QUAL_THRESHOLD', type=float, default=0,
        help='Phred scaled error probability threshold used for genotyping (default: %(default)s). Must be at least 0. '
        'If error probability of genotype is higher, genotype ./. is output.')
    arg('--keep-untagged', action='store_true', 
        help='Consider the untagged reads (reads without a haplotag) in the genotyping. (default: False)')
    

    arg = parser.add_argument_group('Realignment parameters').add_argument
    arg('--realign-mode', metavar='MODE', default="ed",
        help='Select method which will be used to calculate realignment scores. Available methods are: "wfa_full", "wfa_score", and "ed". (default: %(default)s).'
        ' WARNING: "wfa_full" might require a lot of memory if vcf contains large alleles.')
    arg('--overhang', metavar='OVERHANG', default=10, type=int,
        help='When --reference is used, extend alignment by this many bases to left and right when realigning (default: %(default)s).')
    arg('--gap-start', metavar='GAPSTART', default=3, type=float,
        help='gap starting penalty in case wfa is used (default: %(default)s).')
    arg('--gap-extend', metavar='GAPEXTEND', default=1, type=float,
        help='gap extend penalty in case wfa is used (default: %(default)s).')
    arg('--mismatch', metavar='MISMATCH', default=2, type=float,
        help='mismatch cost in case wfa is used (default: %(default)s)')
    
    arg = parser.add_argument_group('Emission probability parameters (Parameters should add up to 1). These are used when mode is "wfa_full"').add_argument
    arg('--match-probability', metavar='MATCH_PROBABILITY', default=0.85, type=float,
        help='probability of match in alignment CIGAR (default: %(default)s)')
    arg('--mismatch-probability', metavar='MISMATCH_PROBABILITY', default=0.05, type=float,
        help='probability of mismatch in alignment CIGAR (default: %(default)s)')
    arg('--insertion-probability', metavar='INSERTION_PROBABILITY', default=0.05, type=float,
        help='probability of insertion in alignment CIGAR (default: %(default)s)')
    arg('--deletion-probability', metavar='DELETION_PROBABILITY', default=0.05, type=float,
        help='probability of deletion in alignment CIGAR (default: %(default)s)')
    
    arg = parser.add_argument_group('HMM parameters').add_argument
    arg('--recombrate', metavar='RECOMBRATE', type=float, default=1.26,
        help='Recombination rate in cM/Mb (used with --ped). If given, a constant recombination '
        'rate is assumed (default: %(default)gcM/Mb).')
    arg('--eff-pop-size', metavar='EFFPOPSIZE', default = 10, type = int,
        help="Parameter for transition probability computing (default: %(default)s)")
# fmt: on


def validate(args, parser):
    if len(args.mapped_read_files) == 0:
        parser.error("No BAM or GAF files found.")
    if args.gt_qual_threshold < 0:
        parser.error("Genotype quality threshold (gt-qual-threshold) must be at least 0.")
    if not args.reference_fasta and not args.rgfa:
        parser.error("No reference found. Please rGFA with GAF files or FASTA with BAM files.")
    if not args.samples:
        parser.error("The sample name has to be provided. Please enter sample name with --sample.")
    if args.match_probability < 0 or args.mismatch_probability < 0 or args.insertion_probability < 0 or args.deletion_probability < 0:
        parser.error("Emission probability parameters cannot be negative.")
    if args.realign_mode not in ["wfa_full", "wfa_score", "ed"]:
        parser.error("Unknown realignment mode detected.")
    

def main(args):
    run_genotype(**vars(args))
