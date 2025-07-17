graph-assemblies/*
    - synthetic assemblies underlying the test GFA

reads-sample<1/2>.fa
    - developer made fasta files with synthetic reads coming from the graph.
    - read names denote the path they originate from.
    - These were made by copy pasting the node sequence. Hence do not contain the small variant information.

reads-sample<1/2>-smallindels-<base error probability>.fa
    - developer made fasta files with synthetic reads coming from the graph and the small indels introduced by the external vcf.
    - read names denote the path and variants they originate from.
    - command used for sample1: `python scripts/generate_reads_from_vcf.py tests/data/truth/prepare-vcf/with-ext.vcf tests/data/fasta/graph-assemblies/REF.fa tests/data/gfa/smallgraph-complete.gfa tests/data/other/sample1-haplotypes.txt 0.001 tests/data/fasta/reads-sample1-smallindels-0.001.fa tests/data/other/sample1-small-var.tsv`
    - command used for sample2: `python scripts/generate_reads_from_vcf.py tests/data/truth/prepare-vcf/with-ext.vcf tests/data/fasta/graph-assemblies/REF.fa tests/data/gfa/smallgraph-complete.gfa tests/data/other/sample2-haplotypes.txt 0.001 tests/data/fasta/reads-sample2-smallindels-0.001.fa tests/data/other/sample2-small-var.tsv`
    - fasta is not reproducible since the code contains random generators.

*.fa.fai
    - index created by `giggles genotype` command for quick access to reads.

reads-reversed.fa
    - set of reverse complement reads made to test gaf processing and reversing reverse complements.