prepare-vcf/
    - other files used for the tests written for `giggles prepare_vcf`

sample<1/2>-haplotypes.txt
    - developer-generated haplotype paths for the samples.
    - they do not contain small variant information added from the external vcf.

sample<1/2>-small-var.tsv
    - the tsv files contains the small variants that the haplotypes passes through.
    - command used for sample1: `python scripts/generate_reads_from_vcf.py tests/data/truth/prepare-vcf/with-ext.vcf tests/data/fasta/graph-assemblies/REF.fa tests/data/gfa/smallgraph-complete.gfa tests/data/other/sample1-haplotypes.txt 0.001 tests/data/fasta/reads-sample1-smallindels-0.001.fa tests/data/other/sample1-small-var.tsv`
    - command used for sample2: `python scripts/generate_reads_from_vcf.py tests/data/truth/prepare-vcf/with-ext.vcf tests/data/fasta/graph-assemblies/REF.fa tests/data/gfa/smallgraph-complete.gfa tests/data/other/sample2-haplotypes.txt 0.001 tests/data/fasta/reads-sample2-smallindels-0.001.fa tests/data/other/sample2-small-var.tsv`
    - these commands are the same for generating the synthetic reads.
    - the script looks for tests/data/other/sample<1/2>-small-var.tsv. if available, it uses that information.
    - if files is not available, it randomly generates the information and writes the file.