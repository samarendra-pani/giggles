# Generating a VCF file with dummy sites list from a FASTA file
# This will be used as the SNP and small indel sites for testing.

import sys
import random

prob_snp = 0.6  # Probability of a site being a SNP
prob_multi_al = 0.2  # Probability of a site having multiple alternate alleles
    
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

# Code to randomize site type between SNP and small indel. Also the number of alternate alleles.
# size is the size of reference allele.
def generate_site(ref_seq, pos, size, is_snp):
    if is_snp:
        assert size == 1, "SNP size must be 1"
        ref = ref_seq[pos]
        if random.random() < prob_multi_al:
            # Generate multiple alternate alleles
            alts = []
            for _ in range(random.randint(2, 3)):
                # Generate a SNP
                alt = random.choice([n for n in 'ACGT' if n != ref and n not in alts])
                alts.append(alt)
            alt = ','.join(alts)
            return ref, alt
        else:
            # Generate a single alternate allele
            alt = random.choice([n for n in 'ACGT' if n != ref])
            return ref, alt
    else:
        # Generate a small indel
        ref = ref_seq[pos:pos+size]
        if random.random() < prob_multi_al:
            # Generate multiple alternate alleles
            alts = []
            for _ in range(random.randint(1, 3)):
                alt_size = random.randint(1, 4)  # Indel size between 1 and 3
                alts.append(''.join(random.choices('ACGT', k=alt_size)))
            alt = ','.join(list(set(alts)))
            return ref, alt
        else:
            alt_size = random.randint(1, 4)  # Indel size between 1 and 3
            if alt_size == 1:
                alt = random.choice([n for n in 'ACGT' if n != ref[0]])
            elif alt_size == 2:
                alt = random.choice([n for n in 'ACGT' if n != ref[0]]) + random.choice([n for n in 'ACGT' if n != ref[-1]])
            else:
                alt = random.choice([n for n in 'ACGT' if n != ref[0]]) + ''.join(random.choices('ACGT', k=alt_size-2)) + random.choice([n for n in 'ACGT' if n != ref[-1]])
            return ref, alt
    
def generate_vcf_sites(fasta, probability):
    print("##fileformat=VCFv4.2")
    print("##source=generate_vcf_sites.py")
    print("##FILTER=<ID=PASS,Description='All filters passed'>")
    print("##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>")
    for seq_id, seq in fasta.items():
        length = len(seq)
        # Print the contig line for each sequence
        print(f"##contig=<ID={seq_id},length={length}>")
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    
    for seq_id, seq in fasta.items():
        length = len(seq)
        num_sites = int(length * probability)  # Number of sites to generate based on the probability
        pos = []
        for _ in range(num_sites):
            pos.append(random.randint(0, length - 1))
        pos = sorted(set(pos))
        
        for p in pos:   
            is_snp = False
            if random.random() < prob_snp:
                size = 1
                is_snp = True
            else:
                if random.random() < 0.5:
                    size = random.randint(2, 4)
                else:
                    size = 1
            ref, alt = generate_site(seq, p, size, is_snp)
            print(f"{seq_id}\t{p+1}\t{'SNP' if is_snp else 'INDEL'}\t{ref}\t{alt}\t.\tPASS\t.")
   
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python generate_vcf_sites.py <fasta_file> <probability>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    probability = float(sys.argv[2])
    
    fasta = read_fasta(fasta_file)
    generate_vcf_sites(fasta, probability)