"""
Genotype variants

Runs only the genotyping algorithm. Genotype Likelihoods are computed using the
forward backward algorithm.
"""

# Code modified from WhatsHap (https://github.com/whatshap/whatshap)

import sys
import platform
import math
# from multiprocessing import Pool
from functools import partial
from contextlib import ExitStack

from giggles import __version__
from giggles.logger import logger
from giggles.vcf import VcfReader, GenotypeVcfWriter
from giggles.core import (
    GenotypeHMM,
    PhredGenotypeLikelihoods,
    Genotype
)
from giggles.utils import (
    UniformRecombinationCostComputer,
)
from giggles.gaf import rGFA
from giggles.timer import StageTimer
from giggles.cli import log_memory_usage
from giggles.utils import select_reads, bin_coeff, determine_genotype
from giggles.cli import ReadSetCreator, read_haplotags


timers = StageTimer()

def genotype_chromosome(variant_table, 
        readset_creator_arguments, 
        haplotags, 
        keep_untagged, 
        max_coverage, 
        gt_prob, 
        recombination_cost_computer, 
        n_haplotypes, 
    ):
    
    chromosome = variant_table.chromosome
    with timers("phased_input_reader"):
        readset_creator = ReadSetCreator(*readset_creator_arguments[0], **readset_creator_arguments[1])

    # create a mapping of genome positions to indices
    var_pos_to_ind = dict()
    n_allele_position = dict()
    allele_references = dict()
    for i in range(len(variant_table.variants)):
        var_pos_to_ind[variant_table.variants[i].position_on_ref] = i
        v = variant_table.variants[i]
        n_allele_position[v.position_on_ref] = len(v.alternative_allele)+1      ##Contains the number of alleles at every variant position
        allele_references[v.position_on_ref] = v.allele_origin        
    
    #Prior genotyping with equal probabilities
    variant_table.query_set_genotype_likelihoods_of(
        [PhredGenotypeLikelihoods([1/(bin_coeff(n_allele_position[pos] + 1, n_allele_position[pos] - 1))] * (bin_coeff(n_allele_position[pos] + 1, n_allele_position[pos] - 1)) , 2, n_allele_position[pos]) for pos in list(var_pos_to_ind.keys())]
    )
    
    # Get the reads belonging to each sample
    with timers("read_alignment"):
        readset = readset_creator.read(
            chromosome, variant_table.variants, haplotags, keep_untagged
        )
    if len(readset) == 0:
        logger.info(f"Skipping chromosome {chromosome} because no reads were found.")
        return

    with timers("select"):
        if max_coverage == None:
            selected_reads = readset
            logger.info(f"Kept {len(readset)} reads in {chromosome}")
        else:
            readset = readset.subset(
                [i for i, read in enumerate(readset) if len(read) >= 2]
            )
            logger.info(f"Kept {len(readset)} reads that cover at least two variants each in {chromosome}")
            selected_reads = select_reads(readset, max_coverage)
    
    # Sorting selected reads
    with timers("alignment_sorting"):
        for read in selected_reads:
            if not read.is_sorted():
                read.sort()
        selected_reads.sort()
        
    # Determine which variants can (in principle) be phased
    accessible_positions = list(var_pos_to_ind.keys())
    accessible_positions_n_allele = []
    accessible_positions_allele_references = []
    for position in accessible_positions:
        accessible_positions_n_allele.append(n_allele_position[position])
        allele_reference_to_list = []
        for ref_sample in allele_references[position]:
            for hap in ref_sample:
                try:
                    allele_reference_to_list.append(int(hap))
                except TypeError:
                    allele_reference_to_list.append(-1)
        accessible_positions_allele_references.append(allele_reference_to_list)
    logger.info(f"Variants in {chromosome} covered by at least one read after read selection: {len(selected_reads.get_positions())}")

    recombination_costs = recombination_cost_computer.compute(accessible_positions)
    
    # Finally, run genotyping algorithm
    with timers("genotyping"):
        forward_backward_table = GenotypeHMM(
            selected_reads,
            recombination_costs,
            n_haplotypes,
            accessible_positions,
            accessible_positions_n_allele,
            accessible_positions_allele_references
        )
        
        # store results
        likelihood_list = variant_table.query_genotype_likelihoods_of()
        genotypes_list = variant_table.query_genotypes_of()

        for pos in range(len(accessible_positions)):
            likelihoods = forward_backward_table.get_genotype_likelihoods(pos, accessible_positions_n_allele[pos])
            # compute genotypes from likelihoods and store information
            geno = determine_genotype(likelihoods, gt_prob, accessible_positions_n_allele[pos])
            assert isinstance(geno, Genotype)
            genotypes_list[var_pos_to_ind[accessible_positions[pos]]] = geno
            likelihood_list[var_pos_to_ind[accessible_positions[pos]]] = likelihoods

        variant_table.query_set_genotypes_of(genotypes_list)
        variant_table.query_set_genotype_likelihoods_of(likelihood_list)

    return (variant_table, chromosome)


def run_genotype(
    alignment_files,
    variant_file,
    rgfa,
    read_fasta_files,
    haplotag_tsv=None,
    keep_untagged=False,
    output=sys.stdout,
    sample="sample",
    # cores=1,
    chromosomes=None,
    mapping_quality=20,
    max_coverage=None,
    gt_qual_threshold=0,
    realign_mode="edit",
    overhang=10,
    gap_start=3,
    gap_extend=1,
    mismatch=2,
    em_reg_constant=10,
    em_score_to_prob_base_constant=math.e,
    match_probability=0.85,
    mismatch_probability=0.05,
    insertion_probability=0.05,
    deletion_probability=0.05,
    recombrate=1.26,
    eff_pop_size = 10
):
    logger.info(f"This is Giggles (genotyping) {__version__} running under Python {platform.python_version()}\n")
    logger.info('======= Working Files')
    logger.info(f"Alignment files: {','.join(alignment_files)}")
    logger.info(f"Read FASTA files: {','.join(read_fasta_files)}")
    logger.info(f"Haplotags for alignments: {'None provided' if haplotag_tsv is None else haplotag_tsv}")
    logger.info(f"Reference GFA file: {rgfa}")
    logger.info(f"Variant sites files: {variant_file}")
    logger.info(f"Writing to: {'Standard output' if output is sys.stdout else output}\n")
    command_line = "(giggles {}) {}".format(__version__, " ".join(sys.argv[1:]))
    with ExitStack() as stack:
        # read the given input files (BAMs, VCFs, ref...)
        if rgfa is not None:
            rgfa = rGFA(reference_path=rgfa)
        readset_creator_args = (alignment_files, rgfa, read_fasta_files)
        readset_creator_kwargs = {'mapq_threshold': mapping_quality,
                'realign_mode': realign_mode,
                'overhang': overhang,
                'gap_start': gap_start,
                'gap_extend': gap_extend,
                'default_mismatch': mismatch,
                'em_prob_params': [match_probability, mismatch_probability, insertion_probability, deletion_probability],
                'reg_const': em_reg_constant,
                'base_const': em_score_to_prob_base_constant}
        
        # reading haplotags
        haplotags = read_haplotags(haplotag_tsv)
        
        # vcf writer for final genotype likelihoods
        vcf_writer = stack.enter_context(GenotypeVcfWriter(command_line=command_line, in_path=variant_file, out_file=output, sample=sample))
        
        # The samples in the gaf or bam is given as input since it will be used to make variant tables with those samples.
        # The variant tables are then simply just updated after the HMM is run.
        vcf_reader = stack.enter_context(
            VcfReader(
                variant_file, indels=True, genotype_likelihoods=False, phases=False, ignore_genotypes=True, required_chr=chromosomes
            )
        )
        recombination_cost_computer = UniformRecombinationCostComputer(recombrate, eff_pop_size)
        # compute genotype likelihood threshold
        gt_prob = 1.0 - (10 ** (-gt_qual_threshold / 10.0))

        # Count number of samples are present in the multisample reference graph vcf file
        # I assume that sample names starting with HG, NA or GM are diploid and others are haploid.
        # TODO: Need a better way to find haploids and diploids. 
        n_haplotypes = 0
        for vcf_reader_sample in list(vcf_reader._vcf_reader.header.samples):
            if vcf_reader_sample[0:2] in ["HG", "NA", "GM"]:
                n_haplotypes += 2
            else:
                n_haplotypes += 1
        
        # creating partial function to give Pool.imap()
        # Note: The phased_input_reader is common to all the processes but python multiprocessing has trouble with pysam-based objects.
        partial_genotype_chromosome = partial(genotype_chromosome, 
            readset_creator_arguments=[readset_creator_args, readset_creator_kwargs], 
            haplotags=haplotags, 
            keep_untagged=keep_untagged, 
            max_coverage=max_coverage, 
            gt_prob=gt_prob, 
            recombination_cost_computer=recombination_cost_computer, 
            n_haplotypes=n_haplotypes)
        
        # No parallel processing
        for variant_table in timers.iterate("parse_vcf", vcf_reader):
            result = partial_genotype_chromosome(variant_table=variant_table)
            if result == None:
                continue
            variant_table, chromosome = result
            # writing VCF outside of parallel processing
            with timers("write_vcf"):
                logger.info(f"======== Writing {chromosome} records")
                vcf_writer.write_genotypes(chromosome, variant_table, indels=True)
        
        '''
        # TODO
        # creating a pool of processes to run the genotyping algorithm in parallel
        # This uses too much RAM (which the memory usage log is also unable to track)
        with Pool(cores) as pool:
            # using pool.imap to yield single values from iterator instead of creating a list with all variant_tables
            results = pool.imap(partial_genotype_chromosome, timers.iterate("parse_vcf", vcf_reader))
            for result in results:
                if result == None:
                    continue
                variant_table, chromosome = result
                # writing VCF outside of parallel processing
                with timers("write_vcf"):
                    logger.info(f"======== Writing {chromosome} records")
                    #vcf_writer.write_genotypes(chromosome, variant_table, indels=True)

                logger.debug(f"Chromosome {chromosome} finished")
        '''

    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info("Time spent reading alignments:               %9.2f s", timers.elapsed("read_alignment"))
    logger.info("Time spent parsing VCF:                      %9.2f s", timers.elapsed("parse_vcf"))
    logger.info("Time spent selecting reads:                  %9.2f s", timers.elapsed("select"))
    logger.info("Time spent sorting selected reads:           %9.2f s", timers.elapsed("alignment_sorting"))
    logger.info("Time spent genotyping:                       %9.2f s", timers.elapsed("genotyping"))
    logger.info("Time spent writing VCF:                      %9.2f s", timers.elapsed("write_vcf"))
    logger.info("Time spent on rest:                          %9.2f s", total_time - timers.sum())
    logger.info("Total elapsed time:                          %9.2f s", total_time)
    logger.info("Total elapsed time:                          %9.2f hr", total_time/3600)


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('variant_file', metavar='VCF', help='VCF file with variants to be genotyped (can be bgzip-compressed).')
    arg('rgfa', metavar='rGFA', help='reference GFA for the GAF alignment (can be bgzip-compressed).')
    arg('alignment_files', metavar='READS', help='Comma separated list of GAF files along with their indexes generated by gaftools sort (can be bgzip-compressed).')
    arg('read_fasta_files', metavar='FASTA', help='Comma separated list of FASTA file with the reads used in GAF files. Please provide in the same order as GAF files. If no index (.fai) exists, it will be created.')
    
    arg('-o', '--output', default=sys.stdout,
        help='Output VCF file. Add .gz to the file name to get compressed output. '
        'If omitted, use standard output.')
    arg('--haplotag-tsv', metavar='HAPLOTAG', 
        help='Comma separated list of TSV file containing the haplotag and phaseset information. Please provide in the same order as GAF files.')
    arg('--sample', dest='sample', metavar='SAMPLE', default='sample',
        help='Name of the sample being genotyped. (default: %(default)s)')
    # arg('--cores', dest='cores', metavar='CORES', default=1, type=int,
    #     help='Number of parallel cores to use for multiprocessing (default: %(default)s)')
    

    arg = parser.add_argument_group('Input pre-processing, selection and filtering').add_argument
    arg('--max-coverage', '-H', metavar='MAX_COV', default=None, type=int,
        help='Reduce coverage to at most MAX_COV. (default: %(default)s)')
    arg('--mapping-quality', '--mapq', metavar='QUAL',
        default=20, type=int, help='Minimum mapping quality (default: %(default)s).')
    arg('--chromosome', dest='chromosomes', metavar='CHROMOSOME', default=[], action='append',
        help='Name of chromosome to genotyped. If not given, all chromosomes in the '
        'input VCF are genotyped. Can be used multiple times.')
    arg('--gt-qual-threshold', metavar='GT_QUAL_THRESHOLD', type=float, default=0,
        help='Phred scaled error probability threshold used for genotyping (default: %(default)s). Must be at least 0. '
        'If error probability of genotype is higher, genotype ./. is output.')
    arg('--keep-untagged', action='store_true', 
        help='Consider the untagged reads (reads without a haplotag) in the genotyping. (default: False)')
    

    arg = parser.add_argument_group('Realignment parameters').add_argument
    arg('--realign-mode', metavar='MODE', default="edit",
        help='Select method which will be used to calculate realignment scores. Available methods are: "wfa_full", "wfa_score", and "edit". (refer to README for more details) (default: %(default)s).')
    arg('--overhang', metavar='OVERHANG', default=10, type=int,
        help='When --reference is used, extend alignment by this many bases to left and right when realigning (default: %(default)s).')
    arg('--gap-start', metavar='GAPSTART', default=3, type=float,
        help='gap starting penalty in case wfa is used (default: %(default)s).')
    arg('--gap-extend', metavar='GAPEXTEND', default=1, type=float,
        help='gap extend penalty in case wfa is used (default: %(default)s).')
    arg('--mismatch', metavar='MISMATCH', default=2, type=float,
        help='mismatch cost in case wfa is used (default: %(default)s)')
    
    arg = parser.add_argument_group('Emission parameters').add_argument
    arg('--em-reg-constant', metavar='REG_CONSTANT', default=10, type=int,
        help='Individual read emissions will be set at a minimum of 10^-REG_CONSTANT (default: %(default)s).')
    arg('--em-score-to-prob-base-constant', metavar='BASE_CONSTANT', default=math.e, type=float,
        help='base value used to conversion of log probabilities to probabilities (default: %(default)s).')
    

    arg = parser.add_argument_group('CIGAR processing parameters (Parameters should add up to 1). These are used when mode is "wfa_full"').add_argument
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
    args.alignment_files = args.alignment_files.split(",")
    args.read_fasta_files = args.read_fasta_files.split(",")
    args.haplotag_tsv = args.haplotag_tsv.split(",") if args.haplotag_tsv else None
    if args.haplotag_tsv is not None and len(args.haplotag_tsv) != len(args.alignment_files):
        parser.error("The number of haplotag TSV files must match the number of GAF files.")
    if not all(f.endswith(".gaf") or f.endswith(".gaf.gz") for f in args.alignment_files):
        parser.error("Only GAF files are supported.")
    if len(args.alignment_files) == 0:
        parser.error("At least one GAF file must be provided.")
    if len(args.alignment_files) != len(args.read_fasta_files):
        parser.error("The number of GAF files must match the number of FASTA files.")
    if args.gt_qual_threshold < 0:
        parser.error("Genotype quality threshold (gt-qual-threshold) must be at least 0.")
    if args.match_probability < 0 or args.mismatch_probability < 0 or args.insertion_probability < 0 or args.deletion_probability < 0:
        parser.error("CIGAR processing parameters cannot be negative.")
    if args.realign_mode not in ["wfa_full", "wfa_score", "edit"]:
        parser.error("Unknown realignment mode detected.")
    
    if args.haplotag_tsv is None:
        if not args.keep_untagged:
            # If no haplotype TSV is provided, we read all alignments
            args.keep_untagged = True
    

def main(args):
    run_genotype(**vars(args))
