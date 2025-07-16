graphaligner-gaf-test.tsv
    - For the test_gaf_coord_on_ref() of test_gaf.py
    - Manually made as truth data

reads-reversed-truth.fa
    - the forward oriented reads corresponding to `tests/data/genotyping/fasta/reads-reversed.fa`

reads-reversed-truth.gaf
    - run using GraphAligner to create the forward oriented GAF alignments
    - GraphAligner version 1.0.20
    - command: `GraphAligner -x vg -g ../gfa/smallgraph-complete.gfa -f reads-reversed-truth.fa -a reads-reversed-truth.gaf`
    - NOTE: GraphAligner will create the alignments but they might not match properly with the reads-reversed.gaf since there are off-by-one errors and incomplete alignments. So requires some manual tweaking.

<file>.sorted.gaf
    - gaftools version 1.2.0
    - command: `gaftools sort --outgaf <file>.sorted.gaf <file>.gaf ../gfa/smallgraph-complete.gfa`
    - also creates <file>.sorted.gaf.gsi