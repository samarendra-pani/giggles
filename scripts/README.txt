extract_seqeuence_around_var.py
    - This is used for pinpointing the variants in the FASTA.
    - Especially useful for small variants since it gives the flanking sequence.

generate_reads_from_vcf.py
    - This is used to create synthetic reads from a VCF and GFA based on the haplotype path your want for the samples.
    - Adjustable parameters for read length and base error probability.

generate_vcf_sites.py
    - Generating a VCF file with dummy sites list from a FASTA file.
    - Adjustable parameters for denisty of variants.

plot_emission_probabilities.py
    - Plotting the trends in emission probabilities varying with temperature.