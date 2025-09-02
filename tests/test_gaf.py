'''
Testing the GAF processing
'''

from giggles.gaf import GafAlignment, rGFA
from giggles.variants import GAFReader
from giggles.vcf import VcfReader
from giggles.cli import PhasedInputReader
from pysam import FastaFile

def parse_output(filename):
    def parse_line(line):
        return tuple(s.strip() for s in line.split("\t"))

    return [parse_line(l) for l in open(filename)]

# testing gaf alignment reversal in case of reverse complement alignment
def test_gaf_reversal(tmp_path):
    
    revcomp_gaf = 'tests/data/gaf/reads-reversed.sorted.gaf'
    revcomp_fa = 'tests/data/fasta/reads-reversed.fa'
    forward_gaf = 'tests/data/truth/reads-reversed-truth.sorted.gaf'
    forward_fa = 'tests/data/truth/reads-reversed-truth.fa'

    graph = 'tests/data/gfa/smallgraph-complete.gfa'

    graph_reader = rGFA(graph)
    revcomp_fasta_reader = FastaFile(revcomp_fa)
    forward_fasta_reader = FastaFile(forward_fa)

    revcomp_alignments = {}
    forward_alignments = {}
    for line in open(revcomp_gaf):
        if line.startswith('#'):
            continue
        alignment = GafAlignment(line=line, source_id=0, fasta=revcomp_fasta_reader)
        alignment, is_reversed = GafAlignment.check_reverse(alignment, graph_reader)
        assert is_reversed
        revcomp_alignments[alignment.read_id] = alignment   
    for line in open(forward_gaf):
        if line.startswith('#'):
            continue
        alignment = GafAlignment(line=line, source_id=0, fasta=forward_fasta_reader)
        alignment, is_reversed = GafAlignment.check_reverse(alignment, graph_reader)
        assert not is_reversed
        forward_alignments[alignment.read_id] = alignment
    
    for read_id in forward_alignments.keys():
        revcomp_alignment = revcomp_alignments[read_id+'_revcomp']
        forward_alignment = forward_alignments[read_id]
        assert revcomp_alignment.q_len == forward_alignment.q_len, (revcomp_alignment.q_len, forward_alignment.q_len)
        assert revcomp_alignment.q_start == forward_alignment.q_start, (revcomp_alignment.q_start, forward_alignment.q_start)
        assert revcomp_alignment.q_end == forward_alignment.q_end, (revcomp_alignment.q_end, forward_alignment.q_end)
        assert revcomp_alignment.p_len == forward_alignment.p_len, (revcomp_alignment.p_len, forward_alignment.p_len)
        assert revcomp_alignment.p_start == forward_alignment.p_start, (revcomp_alignment.p_start, forward_alignment.p_start)
        assert revcomp_alignment.p_end == forward_alignment.p_end, (revcomp_alignment.p_end, forward_alignment.p_end)
        assert revcomp_alignment.path == forward_alignment.path, (revcomp_alignment.path, forward_alignment.path)
        assert revcomp_alignment.cigar == forward_alignment.cigar, (revcomp_alignment.cigar, forward_alignment.cigar)
        

# testing gaf coordinates on reference
def test_gaf_coord_on_ref(tmp_path):
    
    def parse_truth_tsv(file):
        table = {}
        skipped_reads = []
        for line in open(file):
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            if line[0] in table or line[0] in skipped_reads:
                try:
                    table.pop(line[0])
                except KeyError:
                    pass
                continue
            table[line[0]] = [int(line[1]), int(line[2]), int(line[3])]
        return table
    
    graph = 'tests/data/gfa/smallgraph-complete.gfa'
    reads = 'tests/data/fasta/reads-sample1.fa'
    # first testing with graphaligner output
    alignmets = 'tests/data/gaf/reads-sample1-GA.sorted.gaf'
    truth_table = 'tests/data/truth/graphaligner-gaf-test.tsv'
    
    graph_reader = rGFA(graph)
    fasta_reader = FastaFile(reads)

    truth = parse_truth_tsv(truth_table)
    for line in open(alignmets):
        alignment = GafAlignment(line=line, source_id=0, fasta=fasta_reader)
        if 'revcomp' in alignment.read_id:
            continue
        
        alignment, is_reversed = GafAlignment.check_reverse(alignment, graph_reader)

        alignment_start_on_ref = GafAlignment.get_alignment_start_on_ref(alignment=alignment, rgfa=graph_reader)
        alignment_end_on_ref = GafAlignment.get_alignment_end_on_ref(alignment=alignment, rgfa=graph_reader)

        if alignment.read_id not in truth:
            continue

        if is_reversed:
            assert truth[alignment.read_id][0] == 1
        else:
            assert truth[alignment.read_id][0] == 0
        assert alignment_start_on_ref[0] == truth[alignment.read_id][1], (alignment_start_on_ref[0], truth[alignment.read_id][1])
        assert alignment_end_on_ref[0] == truth[alignment.read_id][2], (alignment_end_on_ref[0], truth[alignment.read_id][2])


# testing the identification of variant position and lengths on the gaf alignment paths
def test_gaf_variant_position_identification(tmp_path):
    def parse_truth_json(file):
        import json
        with open(file) as f:
            return json.load(f)

    graph = 'tests/data/gfa/smallgraph-complete.gfa'
    variant_file = 'tests/data/truth/prepare-vcf/with-ext.vcf'
    reads = 'tests/data/fasta/var-position-testing-reads.fa'
    alignments = 'tests/data/gaf/var-position-testing.sorted.gaf'
    truth_file = parse_truth_json('tests/data/truth/var-position-testing.json')

    graph_reader = rGFA(graph)
    gaf_reader = GAFReader(paths=[alignments], reference=graph_reader, read_fasta=reads)
    vcf_reader = VcfReader(variant_file, indels=True, genotype_likelihoods=False, phases=False, ignore_genotypes=True, required_chr=None)
    for variant_table in vcf_reader:
        chromosome = variant_table.chromosome
        alignments = gaf_reader._usable_alignments(chromosome=chromosome)
        updated_variants = gaf_reader._update_variants_in_alignments(alignments=alignments, variants=variant_table.variants)
        for variants_in_alignment, alignment, _, _, partial_bubble_start_of_chrom, partial_bubble_end_of_chrom in updated_variants:
            read_id = alignment.read_id
            start_on_ref, _, start_scaffold = GafAlignment.get_alignment_start_on_ref(alignment=alignment, rgfa=graph_reader)
            end_on_ref, _, end_scaffold = GafAlignment.get_alignment_end_on_ref(alignment=alignment, rgfa=graph_reader)
            if read_id not in truth_file:
                continue
            truth_info = truth_file[read_id]

            assert start_on_ref == truth_info['alignment_start_on_ref'], f"Expected {truth_info['alignment_start_on_ref']} but got {start_on_ref} on read {read_id}"
            assert end_on_ref == truth_info['alignment_end_on_ref'], f"Expected {truth_info['alignment_end_on_ref']} but got {end_on_ref} on read {read_id}"
            assert start_scaffold == truth_info['start_scaffold'], f"Expected {truth_info['start_scaffold']} but got {start_scaffold} on read {read_id}"
            assert end_scaffold == truth_info['end_scaffold'], f"Expected {truth_info['end_scaffold']} but got {end_scaffold} on read {read_id}"
            assert partial_bubble_start_of_chrom == truth_info['partial_bubble_start_of_chrom'], f"Expected {truth_info['partial_bubble_start_of_chrom']} but got {partial_bubble_start_of_chrom} on read {read_id}"
            assert partial_bubble_end_of_chrom == truth_info['partial_bubble_end_of_chrom'], f"Expected {truth_info['partial_bubble_end_of_chrom']} but got {partial_bubble_end_of_chrom} on read {read_id}"
            assert len(variants_in_alignment) == truth_info['number_variants_in_alignment']
            for index, variant in enumerate(variants_in_alignment):
                if variant.id in truth_info['variants_in_alignment']:
                    assert index == truth_info['variants_in_alignment'][variant.id]['index'], f"Expected {truth_info['variants_in_alignment'][variant.id]['index']} but got {index} on variant {variant.id} on read {read_id}"
                    assert variant.position == truth_info['variants_in_alignment'][variant.id]['position'], f"Expected {truth_info['variants_in_alignment'][variant.id]['position']} but got {variant.position} on variant {variant.id} on read {read_id}"
                    assert variant.length_on_path == truth_info['variants_in_alignment'][variant.id]['length_on_path'], f"Expected {truth_info['variants_in_alignment'][variant.id]['length_on_path']} but got {variant.length_on_path} on variant {variant.id} on read {read_id}"


'''
# testing the identification of variants in an alignment
def test_gaf_variant_identification(tmp_path):
    
    graph = 'tests/data/gfa/smallgraph-complete.gfa'
    variant_file = 'tests/data/truth/prepare-vcf/with-ext.vcf'
    
    reads_1 = 'tests/data/fasta/reads-sample1.fa'
    alignmets_1 = 'tests/data/gaf/reads-sample1-GA.sorted.gaf'
    reads_2 = 'tests/data/fasta/reads-sample2.fa'
    alignmets_2 = 'tests/data/gaf/reads-sample2-GA.sorted.gaf'
    
    graph_reader = rGFA(graph)
    gaf_reader = GAFReader(paths=[alignmets_1], reference=graph_reader, read_fasta=reads_1)
    vcf_reader = VcfReader(variant_file, indels=True, genotype_likelihoods=False, phases=False, ignore_genotypes=True, required_chr=None)

    for variant_table in vcf_reader:
        chromosome = variant_table.chromosome
        readset = gaf_reader.read(chromosome=chromosome, variants=variant_table.variants, reference=None)
        print(readset)
        

    assert False

    variant_pointer = 0
    for line in open(alignmets):
        alignment = GafAlignment(line=line, source_id=0, fasta=fasta_reader)
        if 'revcomp' in alignment.read_id:
            continue
        
        alignment, is_reversed = GafAlignment.check_reverse(alignment, graph_reader)

        variants = GafAlignment.(alignment=alignment, rgfa=graph_reader)

        if is_reversed:
            assert len(variants) == 0
        else:
            assert len(variants) > 0
'''