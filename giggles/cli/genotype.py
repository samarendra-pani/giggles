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
    get_max_genotype_ploidy,
    GenotypingAlgorithm,
    Genotype
)
from giggles.utils import (
    UniformRecombinationCostComputer,
)
from giggles.gaf import rGFA
from giggles.timer import StageTimer
from giggles.cli import log_memory_usage
from giggles.utils import update_reads_with_selected, determine_genotype
from giggles.cli import ReadSetCreator, read_haplotags


timers = StageTimer()

def genotype_chromosome(variant_table, 
        readset_creator, 
        haplotags, 
        keep_untagged, 
        max_coverage, 
        gt_prob, 
        recombination_cost_computer, 
        n_haplotypes,
        ploidy,
        temperature
    ):

    chromosome = variant_table.chromosome
    # create a mapping of genome positions to indices
    var_pos_to_ind = dict()
    positions_list = []
    n_allele_list = []
    allele_references_list = []
    is_sv_list = []
    variant_count = 0
    logger.info("Collating variant information to pass into C++ core.")
    for i in range(len(variant_table.variants)):
        v = variant_table.variants[i]
        if v.position_on_ref in var_pos_to_ind:
            raise RuntimeError(f'Position {v.position_on_ref} has multiple variant lines.')
        var_pos_to_ind[v.position_on_ref] = i
        positions_list.append(i)
        n_allele_list.append(len(v.alternative_allele)+1)
        is_sv_list.append(v.is_sv())
        allele_reference_to_list = []
        for ref_sample in v.allele_origin:
            for hap in ref_sample:
                try:
                    allele_reference_to_list.append(int(hap))
                except TypeError:
                    allele_reference_to_list.append(-1)
        allele_references_list.append(allele_reference_to_list)
        variant_count += 1
    
    logger.info("Computing recombination costs.")
    recombination_costs = recombination_cost_computer.compute(positions_list)

    #Prior genotyping with equal probabilities
    variant_table.query_set_genotype_likelihoods_of(
        [None]*variant_count
    )
    
    # Get the reads
    with timers("read_alignment"):
        readset = readset_creator.read(
            chromosome, variant_table.variants, haplotags, keep_untagged
        )
    logger.info(f"Successfully created Readset for {chromosome}. Found {len(readset)} reads covering {len(readset.get_positions())} variants.")
    if len(readset) == 0:
        logger.info(f"Skipping chromosome {chromosome} because no reads were found.")
        return
    
    # Have to do the selection for the phasing algorithm
    with timers("select"):
        #readset = readset.subset(
        #    [i for i, read in enumerate(readset) if len(read) >= 2]
        #)
        #logger.info(f"Kept {len(readset)} reads that cover at least two variants each in {chromosome}")
        logger.info(f"Selecting reads for phasing using maximum coverage of {max_coverage}.")
        update_reads_with_selected(readset, max_coverage)
    
    # Sorting selected reads
    logger.info(f"Sorting the reads")
    for read in readset:
        if not read.is_sorted():
            read.sort()

    # Run genotyping algorithm
    with timers("genotyping-phasing"):
        result = GenotypingAlgorithm(readset, 
            recombination_costs,
            n_haplotypes,
            ploidy,
            temperature,
            positions_list,
            n_allele_list,
            allele_references_list,
            is_sv_list
        )
    
        # store results
        likelihood_list = variant_table.query_genotype_likelihoods_of()
        genotypes_list = variant_table.query_genotypes_of()
        assert len(likelihood_list) == len(positions_list)
        assert len(genotypes_list) == len(positions_list)

        for index, _ in enumerate(positions_list):
            likelihoods = result.get_genotype_likelihoods(index, n_allele_list[index])
            # compute genotypes from likelihoods and store information
            geno = determine_genotype(likelihoods, gt_prob, n_allele_list[index], ploidy)
            assert isinstance(geno, Genotype)
            genotypes_list[index] = geno
            likelihood_list[index] = likelihoods

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
    ploidy=2,
    output=sys.stdout,
    sample="sample",
    # cores=1,
    chromosomes=None,
    mapping_quality=20,
    max_coverage=15,
    gt_qual_threshold=0,
    is_custom_graph=False,
    overhang=10,
    temperature=10.0,
    recombrate=1.26,
    eff_pop_size=10
):
    logger.info(f"This is Giggles (genotyping) {__version__} running under Python {platform.python_version()}.\n")
    logger.info('== Working Files ==')
    logger.info(f"Alignment files: {','.join(alignment_files)}")
    logger.info(f"Read FASTA files: {','.join(read_fasta_files)}")
    logger.info(f"Haplotags for alignments: {'None provided' if haplotag_tsv is None else haplotag_tsv}")
    logger.info(f"Reference GFA file: {rgfa}")
    logger.info(f"Variant sites files: {variant_file}")
    logger.info(f"Writing to: {'Standard output' if output is sys.stdout else output}\n")
    logger.info('== Giggles Genotyping ==')
    command_line = "(giggles {}) {}".format(__version__, " ".join(sys.argv[1:]))
    with ExitStack() as stack:
        # read the given input files (BAMs, VCFs, ref...)
        if rgfa is not None:
            rgfa = rGFA(reference_path=rgfa)
        readset_creator_args = (alignment_files, rgfa, read_fasta_files)
        readset_creator_kwargs = {'mapq_threshold': mapping_quality,
                'is_custom_graph': is_custom_graph,
                'overhang': overhang}
        readset_creator_arguments=[readset_creator_args, readset_creator_kwargs]
        readset_creator = stack.enter_context(ReadSetCreator(*readset_creator_arguments[0], **readset_creator_arguments[1]))

        # reading haplotags
        haplotags = read_haplotags(haplotag_tsv)
        
        # vcf writer for final genotype likelihoods
        logger.info("Initializing VCF writer.")
        vcf_writer = stack.enter_context(GenotypeVcfWriter(command_line=command_line, in_path=variant_file, out_file=output, sample=sample))
        
        # The samples in the gaf or bam is given as input since it will be used to make variant tables with those samples.
        # The variant tables are then simply just updated after the HMM is run.
        logger.info("Initializing VCF reader.")
        vcf_reader = stack.enter_context(
            VcfReader(
                path=variant_file, indels=True, required_chr=chromosomes, is_custom_graph=is_custom_graph
            )
        )
        logger.info("Initializing recombination computer.")
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
        
        partial_genotype_chromosome = partial(genotype_chromosome,
            readset_creator=readset_creator, 
            haplotags=haplotags, 
            keep_untagged=keep_untagged, 
            max_coverage=max_coverage, 
            gt_prob=gt_prob, 
            recombination_cost_computer=recombination_cost_computer, 
            n_haplotypes=n_haplotypes,
            ploidy=ploidy,
            temperature=temperature)
        
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
                logger.info(f"======== Finished writing {chromosome} records")

    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info(f"Time spent reading alignments:               {timers.elapsed('read_alignment'):9.2f} s")
    logger.info(f"Time spent parsing VCF:                      {timers.elapsed('parse_vcf'):9.2f} s")
    logger.info(f"Time spent selecting reads:                  {timers.elapsed('select'):9.2f} s")
    logger.info(f"Time spent genotyping:                       {timers.elapsed('genotyping-phasing'):9.2f} s")
    logger.info(f"Time spent writing VCF:                      {timers.elapsed('write_vcf'):9.2f} s")
    logger.info(f"Time spent on rest:                          {total_time - timers.sum():9.2f} s")
    logger.info(f"Total elapsed time (in seconds):             {total_time:9.2f} s")
    logger.info(f"Total elapsed time (in hours):               {total_time/3600:9.2f} hr")


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('variant_file', metavar='VCF', help='VCF file with variants to be genotyped (can be bgzip-compressed).')
    arg('rgfa', metavar='rGFA', help='reference GFA for the GAF alignment (can be bgzip-compressed).')
    arg('alignment_files', metavar='ALIGNMENTS', help='Comma separated list of GAF files along with their indexes generated by gaftools sort (can be bgzip-compressed).')
    arg('read_fasta_files', metavar='FASTAS', help='Comma separated list of FASTA file with the reads used in GAF files. Please provide in the same order as GAF files. If no index (.fai) exists, it will be created.')
    
    arg('-o', '--output', default=sys.stdout,
        help='Output VCF file. Add .gz to the file name to get compressed output. '
        'If omitted, use standard output.')
    # arg('--rounds', dest='rounds', metavar='ROUNDS', default=2, type=int,
    #     help='Number of phasing-genotyping rounds. (default: %(default)s)')
    arg('--haplotag-tsv', metavar='HAPLOTAG', 
        help='Comma separated list of TSV file containing the haplotag and phaseset information. Please provide in the same order as GAF files.')
    arg('--sample', dest='sample', metavar='SAMPLE', default='sample',
        help='Name of the sample being genotyped. (default: %(default)s)')
    arg('--ploidy', dest='ploidy', metavar='PLOIDY', type=int, default=2,
        help=f'Ploidy of the sample being genotyped. Maximum ploidy supported is {get_max_genotype_ploidy()}. (default: %(default)s)')
    # arg('--cores', dest='cores', metavar='CORES', default=1, type=int,
    #     help='Number of parallel cores to use for multiprocessing (default: %(default)s)')
    

    arg = parser.add_argument_group('Input pre-processing, selection and filtering').add_argument
    arg('--max-coverage', '-H', metavar='MAX_COV', default=15, type=int,
        help='Pre-phasing reads by sub-sampling to MAX_COV. (default: %(default)s)')
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
    arg('--is-custom-graph', action='store_true',
        help='The graph is a custom-made graph where the bubble paths are single nodes corresponding to alleles.')
    #arg('--realignment-bandwidth', metavar='BANDWIDTH', default=30,
    #    help='Set a bandwidth to restrict the realignment process (default: %(default)s).')
    arg('--overhang', metavar='OVERHANG', default=10, type=int,
        help='Extend alignment by this many bases to left and right when realigning (default: %(default)s).')
    arg('--temperature', metavar='TEMPERATURE', default=10.0, type=float,
        help='Parameter to adjust the exponential decay with higher sequence divergence in realignment. High temperature causes '
        'more divergent sequences to decay faster in terms of probability (default: %(default)s).')
    

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
    if args.temperature <= 0:
        parser.error("The temperature parameter as to be positive.")
    if args.haplotag_tsv is not None and len(args.haplotag_tsv) != len(args.alignment_files):
        parser.error("The number of haplotag TSV files must match the number of GAF files.")
    if not all(f.endswith(".gaf") or f.endswith(".gaf.gz") for f in args.alignment_files):
        parser.error("Only GAF files are supported.")
    if len(args.alignment_files) == 0:
        parser.error("At least one GAF file must be provided.")
    if len(args.alignment_files) != len(args.read_fasta_files):
        parser.error("The number of GAF files must match the number of FASTA files.")
    if args.ploidy not in [1, 2]:
        parser.error("Currently only ploidy 1 and 2 are supported.")
    if args.gt_qual_threshold < 0:
        parser.error("Genotype quality threshold (gt-qual-threshold) must be at least 0.")
    
    if args.haplotag_tsv is None:
        if not args.keep_untagged:
            # If no haplotype TSV is provided, we read all alignments
            args.keep_untagged = True
    

def main(args):
    run_genotype(**vars(args))
