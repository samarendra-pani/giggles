# Extract the sequence around a variant in a VCF using a FASTA file.
# For pinpointing the variant in a FASTA files.

import sys

vcf_file = sys.argv[1]
fasta_file = sys.argv[2]
window_size = int(sys.argv[3])
skip_svs = sys.argv[4] == '1'

def read_vcf(vcf_file):
    variants = []
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                chrom = parts[0]
                pos = int(parts[1])
                id = parts[2]
                if skip_svs and '>' in id:
                    continue
                ref = parts[3]
                alts = parts[4].split(',')
                variants.append((chrom, pos, ref, alts))
    return variants

def read_fasta(fasta_file):
    fasta = {}
    with open(fasta_file, 'r') as fa:
        seq_id = None
        for line in fa:
            line = line.strip()
            if line.startswith('>'):
                seq_id = line[1:]  # Remove '>'
                fasta[seq_id] = ''
            else:
                fasta[seq_id] += line
    return fasta

def extract_sequences(variants, fasta):
    for chrom, pos, ref, alts in variants:
        if chrom in fasta:
            sequence = fasta[chrom]
            start = max(0, pos-1 - window_size)  # 10 bases before the variant
            end = min(len(sequence), pos-1 + len(ref) + window_size)  # 10 bases after the variant
            
            print(f"Extracting sequence for variant at {chrom}:{pos} ({ref} -> {alts}). Reference allele in FASTA: {sequence[pos-1:pos-1+len(ref)]}")
            print(f"\tREF: {sequence[start:start+window_size]} | {ref} | {sequence[end-window_size:end]}")
            for alt in alts:
                print(f"\tALT ({alt}): {sequence[start:start+window_size]} | {alt} | {sequence[end-window_size:end]}")

if __name__ == "__main__":
    variants = read_vcf(vcf_file)
    fasta = read_fasta(fasta_file)
    extract_sequences(variants, fasta)