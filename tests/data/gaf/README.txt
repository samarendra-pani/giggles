assembly-alignment/*
    - Alignment of the assemblies found under `tests/data/fasta/graph-assemblies/`
    - GraphAligner version 1.0.20
    - command: `GraphAligner -x vg -g ../../gfa/smallgraph.gfa -f ../../fasta/graph-assemblies/<assembly>.fa -a <assembly>.gaf`

reads-sample1-GA.gaf
    - Alignments of the FASTA `tests/data/fasta/reads-sample1.fa` (does not contain small variant info)
    - GraphAligner version 1.0.20
    - command: `GraphAligner -x vg -g ../gfa/smallgraph-complete.gfa -f ../fasta/reads-sample1.fa -a smallgraph-graphaligner.gaf`

reads-sample<1/2>-smallindels-<base error probability>-GA.gaf
    - Alignments of the FASTA `tests/data/fasta/reads-sample<1/2>-smallindels-<base error probability>.fa`
    - GraphAligner version 1.0.20
    - command for sample1, BEP=0.001: `GraphAligner -x vg -g ../gfa/smallgraph-complete.gfa -f ../fasta/reads-sample1-smallindels-0.001.fa -a reads-sample1-smallindels-0.001-GA.gaf`
    - command for sample2, BEP=0.001: `GraphAligner -x vg -g ../gfa/smallgraph-complete.gfa -f ../fasta/reads-sample2-smallindels-0.001.fa -a reads-sample2-smallindels-0.001-GA.gaf`

reads-reversed.gaf
    - reverse complement reads aligned and checked against alignment called with the same reads but in forward strand.
    - GraphAligner version 1.0.20
    - command: `GraphAligner -x vg -g ../gfa/smallgraph-complete.gfa -f ../fasta/reads-reversed.fa -a reads-reversed.gaf`

<file>.sorted.gaf
    - gaftools version 1.2.0
    - command: `gaftools sort --outgaf <file>.sorted.gaf <file>.gaf ../gfa/smallgraph-complete.gfa`
    - also creates <file>.sorted.gaf.gsi