'''
Testing the GAF processing
'''

from giggles.gaf import GafAlignment, rGFA
from giggles.variants import GAFReader
from giggles.vcf import VcfReader
from giggles.logger import logger
from pysam import FastaFile


logger.set_level('TRACE')

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


# testing the identification of variant position and lengths on the gaf alignment paths.
# tested with the vcf with external variants.
def test_gaf_variant_position_identification_with_ext(tmp_path):
    def parse_truth_json(file):
        import json
        with open(file) as f:
            return json.load(f)

    graph = 'tests/data/gfa/smallgraph-complete.gfa'
    variant_file = 'tests/data/truth/prepare-vcf/with-ext.vcf'
    reads = 'tests/data/fasta/var-position-testing-reads.fa'
    alignments = 'tests/data/gaf/var-position-testing.sorted.gaf'
    truth_data = parse_truth_json('tests/data/truth/var-position-testing-with-ext.json')

    graph_reader = rGFA(graph)
    gaf_reader = GAFReader(alignment_files=[alignments], reference=graph_reader, read_fasta_files=[reads])
    vcf_reader = VcfReader(variant_file, indels=True)
    for variant_table in vcf_reader:
        chromosome = variant_table.chromosome
        alignments = gaf_reader._usable_alignments(chromosome=chromosome)
        updated_variants = gaf_reader._update_variants_in_alignments(alignments=alignments, variants=variant_table.variants)
        for variants_in_alignment, alignment, _ in updated_variants:
            count_not_in_test = 0
            read_id = alignment.read_id
            expected = truth_data[read_id]
            
            # --- TEST BLOCK 1: Alignment Coordinates ---
            # Recalculate these for verification
            calc_start, _, calc_start_node = GafAlignment.get_alignment_start_on_ref(alignment, graph_reader)
            calc_end, _, calc_end_node = GafAlignment.get_alignment_end_on_ref(alignment, graph_reader)

            # Pytest handles the diff, no need for custom messages usually
            assert calc_start == expected['alignment_start_on_ref']
            assert calc_end == expected['alignment_end_on_ref']
            assert calc_start_node == expected['start_scaffold']
            assert calc_end_node == expected['end_scaffold']
            
            if variants_in_alignment is None:
                assert expected['variants_in_alignment'] == None
                assert expected['number_variants_in_alignment'] == 0
                continue

            # --- TEST BLOCK 3: Variant Details ---
            for index, variant in enumerate(variants_in_alignment):
                if variant.id in expected['variants_in_alignment']:
                    expected_var_info = expected['variants_in_alignment'][variant.id]
                else:
                    count_not_in_test += 1
                    continue

                # Consolidate assertions for cleaner failure reports
                # Create a dict of what we found to compare against expected dict
                actual_var_info = {
                    'index': index,
                    'position': variant.position,
                    'length_on_path': variant.length_on_path,
                    'state': variant.state
                }
                
                # Compare the dictionary subset. 
                # This gives a beautiful diff in the terminal if it fails.
                assert actual_var_info == {
                    k: expected_var_info[k] for k in actual_var_info
                }, f"Mismatch in variant details for {variant.id} in read {read_id}"
            
            assert count_not_in_test + len(expected['variants_in_alignment']) == expected['number_variants_in_alignment']
            assert len(variants_in_alignment) == expected['number_variants_in_alignment']



# testing the identification of variant position and lengths on the gaf alignment paths.
# tested with the vcf without external variants.
def test_gaf_variant_position_identification_no_ext(tmp_path):
    def parse_truth_json(file):
        import json
        with open(file) as f:
            return json.load(f)

    graph = 'tests/data/gfa/smallgraph-complete.gfa'
    variant_file = 'tests/data/truth/prepare-vcf/no-ext.vcf'
    reads = 'tests/data/fasta/var-position-testing-reads.fa'
    alignments = 'tests/data/gaf/var-position-testing.sorted.gaf'
    truth_data = parse_truth_json('tests/data/truth/var-position-testing-no-ext.json')

    graph_reader = rGFA(graph)
    gaf_reader = GAFReader(alignment_files=[alignments], reference=graph_reader, read_fasta_files=[reads])
    vcf_reader = VcfReader(variant_file, indels=True)
    for variant_table in vcf_reader:
        chromosome = variant_table.chromosome
        alignments = gaf_reader._usable_alignments(chromosome=chromosome)
        updated_variants = gaf_reader._update_variants_in_alignments(alignments=alignments, variants=variant_table.variants)
        for variants_in_alignment, alignment, _ in updated_variants:
            count_not_in_test = 0
            read_id = alignment.read_id
            expected = truth_data[read_id]
            
            # --- TEST BLOCK 1: Alignment Coordinates ---
            # Recalculate these for verification
            calc_start, _, calc_start_node = GafAlignment.get_alignment_start_on_ref(alignment, graph_reader)
            calc_end, _, calc_end_node = GafAlignment.get_alignment_end_on_ref(alignment, graph_reader)

            # Pytest handles the diff, no need for custom messages usually
            assert calc_start == expected['alignment_start_on_ref']
            assert calc_end == expected['alignment_end_on_ref']
            assert calc_start_node == expected['start_scaffold']
            assert calc_end_node == expected['end_scaffold']
            
            if variants_in_alignment is None:
                assert expected['variants_in_alignment'] == None
                assert expected['number_variants_in_alignment'] == 0
                continue

            # --- TEST BLOCK 3: Variant Details ---
            for index, variant in enumerate(variants_in_alignment):
                if variant.id in expected['variants_in_alignment']:
                    expected_var_info = expected['variants_in_alignment'][variant.id]
                else:
                    count_not_in_test += 1
                    continue

                # Consolidate assertions for cleaner failure reports
                # Create a dict of what we found to compare against expected dict
                actual_var_info = {
                    'index': index,
                    'position': variant.position,
                    'length_on_path': variant.length_on_path,
                    'state': variant.state
                }
                
                # Compare the dictionary subset. 
                # This gives a beautiful diff in the terminal if it fails.
                assert actual_var_info == {
                    k: expected_var_info[k] for k in actual_var_info
                }, f"Mismatch in variant details for {variant.id} in read {read_id}"
            
            assert count_not_in_test + len(expected['variants_in_alignment']) == expected['number_variants_in_alignment']
            assert len(variants_in_alignment) == expected['number_variants_in_alignment']