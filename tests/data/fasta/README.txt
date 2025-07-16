graph-assemblies/*
    - synthetic assemblies underlying the test GFA

reads.fa
    - developer made fasta files with synthetic reads coming from the graph.
    - read names denote the path they originate from.

reads.fa.fai
    - index created by `giggles genotype` command for quick access to reads.

reads-reversed.fa
    - set of reverse complement reads made to test gaf processing and reversing reverse complements.