assembly-alignment/*
    - Alignment of the assemblies found under `tests/data/fasta/graph-assemblies/`
    - GraphAligner version 1.0.20
    - command: `GraphAligner -x vg -g ../../gfa/smallgraph.gfa -f ../../fasta/graph-assemblies/<assembly>.fa -a <assembly>.gaf`

smallgraph-graphaligner.gaf
    - GraphAligner version 1.0.20
    - command: `GraphAligner -x vg -g ../gfa/smallgraph-complete.gfa -f ../fasta/reads.fa -a smallgraph-graphaligner.gaf`

smallgraph-minigraph.gaf
    - minigraph version 0.21
    - command: `minigraph --vc -c -x lr -o smallgraph-minigraph.gaf ../gfa/smallgraph-complete.gfa ../fasta/reads.fa 2> smallgraph-minigraph.log`

reads-reversed.gaf
    - reverse complement reads aligned and checked against alignment called with the same reads but in forward strand.
    - GraphAligner version 1.0.20
    - command: `GraphAligner -x vg -g ../gfa/smallgraph-complete.gfa -f ../fasta/reads-reversed.fa -a reads-reversed.gaf`

<file>.sorted.gaf
    - gaftools version 1.2.0
    - command: `gaftools sort --outgaf <file>.sorted.gaf <file>.gaf ../gfa/smallgraph-complete.gfa`
    - also creates <file>.sorted.gaf.gsi