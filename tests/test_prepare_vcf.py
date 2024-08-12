'''
Testing the prepare-vcf subcommand of giggles
'''

from giggles.cli.prepare_vcf import run

def parse_output(filename):
    def parse_line(line):
        return tuple(s.strip() for s in line.split("\t"))

    return [parse_line(l) for l in open(filename)]

# testing vcf creation with haploid assembly gaf
def test_haploid(tmp_path):
    
    graph = 'tests/data/smallgraph-ordered.gfa'
    hap_fofn = 'tests/data/fofn-1.txt'
    test_out = str(tmp_path)+'/prepared-vcf.vcf'

    true_out = 'tests/data/smallgraph-prepared-vcf-haploid.vcf'

    run(gfa=graph, haploid=hap_fofn, output=test_out)

    test_output_lines = parse_output(test_out)
    true_output_lines = parse_output(true_out)

    assert len(test_output_lines)==len(true_output_lines)
    for n in range(len(test_output_lines)):
        assert test_output_lines[n]==true_output_lines[n]
