'''
Testing the GAF processing
'''

from giggles.gaf import GafAlignment, rGFA
from pysam import FastaFile

def parse_output(filename):
    def parse_line(line):
        return tuple(s.strip() for s in line.split("\t"))

    return [parse_line(l) for l in open(filename)]

# testing gaf alignment reversal in case of reverse complement alignment
def test_gaf_reversal(tmp_path):
    
    revcomp_gaf = 'tests/data/genotyping/gaf/reads-reversed.sorted.gaf'
    revcomp_fa = 'tests/data/genotyping/fasta/reads-reversed.fa'
    forward_gaf = 'tests/data/genotyping/truth/reads-reversed-truth.sorted.gaf'
    forward_fa = 'tests/data/genotyping/truth/reads-reversed-truth.fa'

    graph = 'tests/data/genotyping/gfa/smallgraph-complete.gfa'

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
    
    graph = 'tests/data/genotyping/gfa/smallgraph-complete.gfa'
    reads = 'tests/data/genotyping/fasta/reads.fa'
    # first testing with graphaligner output
    alignmets = 'tests/data/genotyping/gaf/smallgraph-graphaligner.sorted.gaf'
    truth_table = 'tests/data/genotyping/truth/graphaligner-gaf-test.tsv'
    
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