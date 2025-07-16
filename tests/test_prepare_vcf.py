'''
Testing the prepare-vcf subcommand of giggles
'''

from giggles.cli.prepare_vcf import run

def parse_output(filename):
    def parse_line(line):
        return tuple(s.strip() for s in line.split("\t"))

    return [parse_line(l) for l in open(filename)]

# testing vcf creation (without any external small variants)
def test_vcf_creation_no_ext(tmp_path):
    
    graph = 'tests/data/gfa/smallgraph-complete.gfa'
    hap_fofn = 'tests/data/other/prepare-vcf/fofn-haploid.txt'
    dip_fofn = 'tests/data/other/prepare-vcf/fofn-diploid.txt'
    test_out = str(tmp_path)+'/prepared-vcf.vcf'

    true_out = 'tests/data/truth/prepare-vcf/no-ext.vcf'

    run(gfa=graph, haploid=hap_fofn, diploid=dip_fofn, output=test_out)

    test_output_lines = parse_output(test_out)
    true_output_lines = parse_output(true_out)

    assert len(test_output_lines)==len(true_output_lines)
    for n in range(len(test_output_lines)):
        assert test_output_lines[n]==true_output_lines[n]
'''
# testing vcf creation (with external small variants)
def test_vcf_creation_with_ext(tmp_path):
    
    graph = 'tests/data/gfa/smallgraph-complete.gfa'
    hap_fofn = 'tests/data/other/prepare-vcf/fofn-haploid.txt'
    dip_fofn = 'tests/data/other/prepare-vcf/fofn-diploid.txt'
    ext_vcf = 'tests/data/vcf/small-indel.vcf'
    test_out = str(tmp_path)+'/prepared-vcf.vcf'

    true_out = 'tests/data/truth/prepare-vcf/with-ext.vcf'

    run(gfa=graph, haploid=hap_fofn, diploid=dip_fofn, external_vcf=ext_vcf, output=test_out)

    test_output_lines = parse_output(test_out)
    true_output_lines = parse_output(true_out)

    assert len(test_output_lines)==len(true_output_lines)
    for n in range(len(test_output_lines)):
        assert test_output_lines[n]==true_output_lines[n]
'''
# testing correct identification of variants on scaffold nodes
def test_scaffold_variants(tmp_path):

    graph = 'tests/data/gfa/smallgraph-complete.gfa'
    ext_vcf = 'tests/data/vcf/small-indel.vcf'
    truth_out = 'tests/data/truth/prepare-vcf/scaffold-variants.tsv'

    def parse_scaffold_variants(file):
        scaffold_variants = {}
        with open(file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                scaffold_variants[(parts[0], int(parts[1]))] = parts[2] == 'TRUE'
        return scaffold_variants

    from giggles.cli.prepare_vcf import read_gfa, read_external_vcf, label_ext_variants
    from collections import defaultdict, namedtuple

    Node = namedtuple("Node", ["SN", "SO", "SR", "LN", "BO", "NO", "Seq"], defaults=[-1, -1, None, -1, -1, -1, None])
    Edge = namedtuple("Edge", ["next", "prev"], defaults=[None, None])
    nodes = defaultdict(lambda: Node())
    edges = defaultdict(lambda: Edge())
    read_gfa(graph, nodes, edges)
    ext_variants = read_external_vcf(ext_vcf)
    scaffold_nodes = [node_id for node_id, value in sorted(nodes.items(), key=lambda item: item[1].BO) if value.NO==0]
    
    # find the variants that fall on the scaffold nodes
    label_ext_variants(ext_variants, scaffold_nodes, nodes)
    scaffold_variants_truth = parse_scaffold_variants(truth_out)
    for variant in ext_variants['chr1']:
        assert scaffold_variants_truth[('chr1', variant['POS'])] == variant['on_scaffold'], f"Variant at {variant['POS']} on scaffold mismatch: expected {scaffold_variants_truth[('chr1', variant['POS'])]}, got {variant['on_scaffold']}"