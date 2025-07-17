# The script creates reads from a VCF file and generates a FASTA file with the reads.
# Uses the following files
# VCF: vcf file with all the variants
# Reference FASTA: reference genome in FASTA format
# Haplotype list: file with haplotype node path (probably only containing the bubbles. The small variants can be decided by the script)
# Base error probability: probability of a base being incorrect
# FASTA output: output file for the generated reads in FASTA format
# Alternate Allele file: file with alternate alleles. If the file does not exist, it will be written based on the random generation of variants.
# # Average read length: average read length for the reads. If not provided, it will be set to 2000.

import sys
import random
from collections import defaultdict
import re
import logging
import os
import numpy

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)

def read_vcf(vcf_file):
    variants = {}
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                chrom = parts[0]
                pos = int(parts[1])
                id = parts[2]
                ref = parts[3]
                alts = parts[4].split(',')
                alt_dict = {}
                if 'EXT' not in id:
                    # bubble variants
                    # process the INFO field to get the allele traversels.
                    info = parts[7].split(';')
                    for info_part in info:
                        if info_part.startswith('AT='):
                            ats = info_part.split('=')[1].split(',')
                            assert len(ats) == len(alts)+1, "Number of allele traversals should be reference allele + number of alternate alleles"
                            for index, at in enumerate(ats):
                                if index == 0:
                                    # reference allele
                                    alt_dict[at] = ref
                                    continue
                                alt_dict[at] = alts[index-1]
                if chrom not in variants:
                    variants[chrom] = []
                variants[chrom].append((pos, id, ref, alts, alt_dict))
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


def read_gfa(gfa_file):
    nodes = {}
    with open(gfa_file, 'r') as gfa:
        for line in gfa:
            if line.startswith('S'):
                parts = line.strip().split('\t')
                node_id = parts[1]
                seq = parts[2]
                length = len(seq)
                so = None
                sr = None
                bo = None
                no = None
                # read tags
                for tag in parts[3:]:
                    if tag.startswith('SO:'):
                        so = int(tag.split(':')[2])
                    elif tag.startswith('SR:'):
                        sr = int(tag.split(':')[2])
                    elif tag.startswith('BO:'):
                        bo = int(tag.split(':')[2])
                    elif tag.startswith('NO:'):
                        no = int(tag.split(':')[2])
                nodes[node_id] = [seq, length, so, sr, bo, no]
    return nodes


def read_haplotype_list(haplotype_file):
    haplotypes = []
    with open(haplotype_file, 'r') as hf:
        for line in hf:
            haplotypes.append(line.strip())
    return haplotypes


def introduce_noise(reads, base_error_prob):
    noisy_reads = defaultdict(list)
    for read_id, read_seq in reads.items():
        noisy_seq = []
        for base in read_seq:
            if random.random() < base_error_prob:
                # Introduce a random base error
                noisy_base = random.choice(['A', 'T', 'C', 'G'])
                noisy_seq.append(noisy_base)
            else:
                noisy_seq.append(base)
        noisy_reads[read_id] = ''.join(noisy_seq)
    return noisy_reads


def read_alt_alleles(alt_alleles_file):
    alt_alleles = {}
    with open(alt_alleles_file, 'r') as af:
        for line in af:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            id = parts[0]
            alleles = parts[1].split(',')
            hap1_index = int(parts[2])
            hap2_index = int(parts[3])
            alt_alleles[id] = [alleles, hap1_index, hap2_index]
    return alt_alleles


def write_alt_alleles(alt_alleles, alt_alleles_file):
    hap1 = alt_alleles[0]
    hap2 = alt_alleles[1]
    with open(alt_alleles_file, 'w') as af:
        af.write("#id\talleles\thap1\thap2\n")  # Header for the alternate alleles file
        for id, (hap1_alleles, hap1_allele, hap1_index) in hap1.items():
            assert id in hap2, f"ID {id} not found in hap2 alleles."
            hap2_alleles, hap2_allele, hap2_index = hap2[id]
            assert hap1_alleles == hap2_alleles, f"Alleles for ID {id} do not match between haplotypes."
            alleles = hap1_alleles
            af.write(f"{id}\t{alleles}\t{hap1_index}\t{hap2_index}\n")
    logging.info(f"Alternate alleles written to {alt_alleles_file}")


def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))


def generate_sequence_from_variants(variants, gfa_nodes, haplotypes, reference, alternate_alleles):
    variant_pointer = 0
    orient = None
    seq = [{}, {}]  # keeping track of the sequences for the external variants
    alternate_variants_selected = [{}, {}] # keeping track of the alternate alleles for the external variants
    # Iterate through each haplotype
    for hap_index, hap in enumerate(haplotypes):
        logging.info(f"Processing haplotype {hap_index+1}: {hap}")
        # splitting the path based on > and <
        for node in list(filter(None, re.split('(>)|(<)', hap))):
            if node in ['>', '<']:
                orient = node
                continue
            if gfa_nodes.get(node) is None:
                sys.exit(f"Node {node} not found in GFA nodes.")
            logging.info(f"Processing node {node} with orientation {orient}")
            node_seq = gfa_nodes[node][0]
            length = gfa_nodes[node][1]
            so = gfa_nodes[node][2]
            no = gfa_nodes[node][5]
            new_start = so
            if no != 0:
                logging.info(f"Node {node} is not a scaffold node. Adding sequence directly.")
                seq[hap_index][node] = node_seq if orient == '>' else reverse_complement(node_seq)
                continue
            logging.info(f"Node {node} is a scaffold node. Processing external variants in it.")
            end = so + length
            while variant_pointer < len(variants):
                pos, id, ref, alts, _ = variants[variant_pointer]
                if 'EXT' not in id:
                    logging.info(f'\t[1] EXT not in id, skipping variant {id}')
                    # Only external variants
                    variant_pointer += 1
                    logging.info(f'\t[1] Variant pointer incremented to {variant_pointer}')
                    continue
                if pos < so:
                    # Variant is before the current node
                    sys.exit(f"Found variant {id} at position {pos} before the start of node {node} at position {so}.")
                    logging.info(f"\t[2] Variant {id} at position {pos} before the start of node {node} at position {so}.")
                    variant_pointer += 1
                    logging.info(f'\t[2] Variant pointer incremented to {variant_pointer}')
                    continue
                if pos < end:
                    # writing the reference sequence before the external variant start
                    seq[hap_index][(node, id, 'ref', new_start)] = node_seq[new_start-so:pos-1-so]
                    assert node_seq[new_start-so:pos-1-so] == reference['chr1'][new_start: pos-1], f"Reference sequence mismatch at {node} for variant {id}. Expected {reference[new_start:pos-1]}, got {node_seq[new_start-so:pos-1-so]}"
                    if alternate_alleles:
                        # extracting the alt allele info from files
                        alleles = alternate_alleles[id][0]
                        allele_index = int(alternate_alleles[id][1+hap_index])
                        allele = alleles[allele_index]
                    else:
                        # decide which allele to use
                        alleles = [ref]+alts
                        # randomly choose an allele
                        allele_index = random.randint(0, len(alleles)-1)
                        allele = alleles[allele_index]
                        alternate_variants_selected[hap_index][id] = [','.join(alleles), allele, allele_index]
                    seq[hap_index][(node, id, 'allele')] = allele
                    new_start = pos - 1 + len(ref)
                    logging.info(f"\t[3] Variant {id} at position {pos} in node {node}, using allele {allele}.")
                    variant_pointer += 1
                    logging.info(f'\t[3] Variant pointer incremented to {variant_pointer}')
                    continue
                if pos >= end:
                    # Variant is after the current node. Breaking the loop to move to the next node.
                    logging.info(f"\t[4] Variant {id} at position {pos} after the end of node {node} at position {end}. Breaking to next node.")
                    break
                variant_pointer += 1
                logging.info(f'\t[5] Variant pointer incremented to {variant_pointer}')
            seq[hap_index][(node, id, 'ref', new_start)] = reference['chr1'][new_start:end]
        variant_pointer = 0  # Reset variant pointer for the next haplotype
    
    return seq, alternate_variants_selected
        

def find_variants_covered_by_read(start, end, coordinate_table):
    """
    Find variants that are covered by a read based on its start and end coordinates.
    """
    variants_covered = []
    for key, (start_coord, end_coord, *rest) in coordinate_table.items():
        if start > end_coord or end < start_coord:
            continue
        # The read overlaps with the variant
        variants_covered.append((key, start_coord, end_coord))
    # sort the variants by their start coordinate
    variants_covered.sort(key=lambda x: x[1])
    return variants_covered


def create_read_id_from_variants(hap_index, read_index, variants_covered):
    """
    Create a unique read ID based on the haplotype index, read index, and variants covered by the read.
    """
    variant_ids = []
    for variant in variants_covered:
        variant = variant[0]
        variant_id = variant[0] if isinstance(variant, tuple) else variant
        variant_ids.append(variant_id)
    variant_ids_str = '_'.join(str(v) for v in variant_ids)
    return f"hap{hap_index+1}_read{read_index}_{variant_ids_str}"

def create_reads_from_haplotypes(haplotypes, avg_read_length, coordinate_table):
    # create a read length distribution based on the average read length
    reads = {}
    for hap_index, hap in enumerate(haplotypes):
        num_reads = random.randint(40, 50)  # Randomly decide the number of reads to generate
        read_length_distribution = numpy.random.normal(avg_read_length, avg_read_length * 0.3, num_reads)
        # logging.info(f"Generating reads for haplotype {hap_index+1}")
        for read_index in range(num_reads):
            read_length = int(read_length_distribution[read_index])
            start = random.randint(0, len(hap) - read_length)
            end = start + read_length
            read_seq = hap[start:end]
            if read_seq == '':
                logging.warning(f"Generated empty read for haplotype {hap_index+1}, read index {read_index}. Skipping.")
                continue
            variants_covered = find_variants_covered_by_read(start, end, coordinate_table[hap_index])
            read_id = create_read_id_from_variants(hap_index, read_index+1, variants_covered)
            reads[read_id] = read_seq
    return reads
            

if __name__ == "__main__":
    if len(sys.argv) < 8:
        print("Usage: python generate_vcf_sites.py <vcf> <reference fasta> <gfa> <haplotype list> <base error probability> <fasta output> <alternate alleles file> [<average read length>]")
        sys.exit(1)
    
    vcf_file = sys.argv[1]
    fasta_file = sys.argv[2]
    gfa = sys.argv[3]
    haplotype_list_file = sys.argv[4]
    base_error_prob = float(sys.argv[5])
    fasta_output = sys.argv[6]  # Output file for the generated reads in FASTA format
    alternate_alleles_file = sys.argv[7] # File with alternate alleles. If file does not exist, it will be written
    avg_read_length = int(sys.argv[8]) if len(sys.argv) > 8 else 2000  # Default read length if not provided
    
    if not (0 <= base_error_prob <= 1):
        print("Base error probability must be between 0 and 1.")
        sys.exit(1)
    
    # Read the VCF file, reference FASTA, GFA, and haplotype list
    variants = read_vcf(vcf_file)
    reference = read_fasta(fasta_file)
    gfa_nodes = read_gfa(gfa)
    haplotypes = read_haplotype_list(haplotype_list_file)
    # checking if the alternate alleles file exists
    alternate_alleles = read_alt_alleles(alternate_alleles_file) if os.path.exists(alternate_alleles_file) else None
    
    # generating read ids from haplotypes and which variants are in the reads
    seq, alt_dict = generate_sequence_from_variants(variants['chr1'], gfa_nodes, haplotypes, reference, alternate_alleles)
    
    # writing the alternate alleles to the file if it does not exist
    if not os.path.exists(alternate_alleles_file):
        write_alt_alleles(alt_dict, alternate_alleles_file)
    
    # create table of start and end coordinate of variants to keep track of variants covered by each read
    haplotypes = ['', '']
    coordinate_table = [{}, {}] # keeping track of the coordinates of the variants in haplotypes

    for hap_index, seq_dict in enumerate(seq):
        for key, value in seq_dict.items():
            if isinstance(key, tuple):
                # inside a scaffold node
                if 'allele' in key:
                    # The string is a variant allele
                    node, id, _ = key
                    start = len(haplotypes[hap_index])
                    end = start + len(value)
                    haplotypes[hap_index] += value
                    coordinate_table[hap_index][id] = (start, end, node, value)
                else:
                    # The string is a reference sequence
                    node, id, _, start_on_ref = key
                    start = len(haplotypes[hap_index])
                    end = start + len(value)
                    haplotypes[hap_index] += value
                    coordinate_table[hap_index][(node, start_on_ref)] = (start, end)
            else:
                # a non-scaffolf node
                start = len(haplotypes[hap_index])
                end = start + len(value)
                haplotypes[hap_index] += value
                coordinate_table[hap_index][key] = (start, end)
    
    # generating reads from haplotypes
    reads = create_reads_from_haplotypes(haplotypes, avg_read_length, coordinate_table)
   
    # introduce base errors in the reads
    noisy_reads = introduce_noise(reads, base_error_prob)
    # shuffle the reads to simulate random sequencing
    noisy_reads = dict(random.sample(list(noisy_reads.items()), len(noisy_reads)))  # Shuffle the reads
    # write fasta
    fasta_writer = open(fasta_output, 'w')
    for read_id, read_seq in noisy_reads.items():
        print(len(read_seq))
        print(f">{read_id}\n{read_seq}", file=fasta_writer)
    fasta_writer.close()
