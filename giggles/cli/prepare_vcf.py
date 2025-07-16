'''
Create VCFs using assembly-to-graph alignments.
'''

import sys
import logging
from collections import defaultdict, namedtuple, abc
import gzip
import re
from copy import deepcopy

logger = logging.getLogger(__name__)

VariantRecord = namedtuple('VariantRecord', ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])

class VCFWriter:
    def __init__(self, path):
        if path == None:
            self.writer = sys.stdout
        else:
            self.writer = open(path, "w")
    
    def write(self, line):
        line += '\n'
        self.writer.write(line)
        
    def close(self):
        self.writer.close()


class Genotype:
    def __init__(self):
        self.maternal = None
        self.paternal = None
    
    def to_string(self):
        return '%s|%s'%(self.paternal, self.maternal)

class Haplotype:
    def __init__(self):
        self.allele = None
        self._remove = False
    
    def to_string(self):
        return '%s'%(self.allele)

def run(
    gfa=None,
    haploid=None,
    diploid=None,
    keep_all_records=False,
    external_vcf=None,
    output=None
):
    
    Node = namedtuple("Node", ["SN", "SO", "SR", "LN", "BO", "NO", "Seq"], defaults=[-1, -1, None, -1, -1, -1, None])
    Edge = namedtuple("Edge", ["next", "prev"], defaults=[None, None])
    nodes = defaultdict(lambda: Node())
    edges = defaultdict(lambda: Edge())
    read_gfa(gfa, nodes, edges)
    ext_variants = read_external_vcf(external_vcf) if external_vcf else {}
    scaffold_nodes = [node_id for node_id, value in sorted(nodes.items(), key=lambda item: item[1].BO) if value.NO==0]
    
    # find the variants that fall on the scaffold nodes
    label_ext_variants(ext_variants, scaffold_nodes, nodes)
    
    diploid_variants, haploid_variants = None, None
    diploid_haplotype_list, haploid_haplotype_list = [], []
    diploid_contigs, haploid_contigs = [], []

    if diploid:
        diploid_variants, diploid_haplotype_list, diploid_contigs = read_assemblies(diploid, scaffold_nodes, nodes)
    if haploid:
        haploid_variants, haploid_haplotype_list, haploid_contigs = read_assemblies(haploid, scaffold_nodes, nodes)

    if diploid:
        variants = diploid_variants
        if haploid:
            for bub in diploid_variants.keys():
                variants[bub].update(haploid_variants[bub])
    else:
        variants = haploid_variants
    
    diploid_haplotype_list.sort()
    haploid_haplotype_list.sort()
    haplotype_list = diploid_haplotype_list+haploid_haplotype_list
    contigs = list(set(diploid_contigs+haploid_contigs))
    contigs.sort()
    ref_alleles = get_reference_alleles(edges, list(variants.keys()))
    
    writer = VCFWriter(output)
    write_vcf(writer, variants, ref_alleles, haplotype_list, contigs, nodes, keep_all_records, ext_variants)
    writer.close()
    

def label_ext_variants(ext_variants, scaffold_nodes, nodes):
    """
    Label external variants based on their position on the scaffold nodes.
    """

    # extracting the coordinates of scaffold nodes
    scaffold_coordiantes = {}
    for node in scaffold_nodes:
        chrom = nodes[node].SN
        start = nodes[node].SO
        end = start + nodes[node].LN
        if chrom not in scaffold_coordiantes:
            scaffold_coordiantes[chrom] = []
        scaffold_coordiantes[chrom].append((start, end))
    
    # sorting the coordinates for each chromosome
    for chrom in scaffold_coordiantes:
        scaffold_coordiantes[chrom].sort(key=lambda x: x[0])
    
    # assuming the VCF is sorted
    # iterating through the external variants and checking if they fall on the scaffold nodes
    for chrom in ext_variants:
        if chrom not in scaffold_coordiantes:
            logger.error(f"Chromosome {chrom} in external variants not found in scaffold nodes. Skipping.")
            continue
        variants = ext_variants[chrom]
        scaffold_pointer = 0
        for variant in variants:
            pos = variant['POS']
            # Move the scaffold pointer to the right position
            while scaffold_pointer < len(scaffold_coordiantes[chrom]) and pos > scaffold_coordiantes[chrom][scaffold_pointer][1]:
                scaffold_pointer += 1
            # Check if the variant falls on a scaffold node
            if scaffold_pointer < len(scaffold_coordiantes[chrom]) and pos >= scaffold_coordiantes[chrom][scaffold_pointer][0] and pos <= scaffold_coordiantes[chrom][scaffold_pointer][1]:
                variant['on_scaffold'] = True
            else:
                variant['on_scaffold'] = False   


def read_external_vcf(external_vcf):
    """
    Read an external VCF file and return a dictionary of variants.
    """
    variants = {}
    with open(external_vcf, 'r') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alts = fields[4].split(',')
            if chrom not in variants:
                variants[chrom] = []
            variants[chrom].append({
                'POS': pos,
                'REF': ref,
                'ALT': alts,
                'on_scaffold': None
            })
    return variants

def read_gfa(gfa, node, edges):
    
    logger.info("\nParsing GFA file and reading sort key information")
    logger.info("\tRequires the following tags: SN, SO, SR, BO, and NO. Either the sequence or the LN tag has to be provided also.")
    if is_file_gzipped(gfa):
        reader = gzip.open(gfa, 'rt')
    else:
        reader = open(gfa, 'r')
    total_nodes = 0
    tagged_nodes = 0
    while True:
        line = reader.readline()
        if not line:
            break
        if line[0] != 'S':
            continue
        total_nodes += 1
        fields = line.split("\t")
        LN = -1
        SN = None
        SO = -1
        SR = -1
        BO = -1
        NO = -1
        for f in fields:
            if f.startswith("LN:i:"):
                LN = int(f[5:])
            if f.startswith("SN:Z:"):
                SN = f[5:]
            if f.startswith("SO:i:"):
                SO = int(f[5:])
            if f.startswith("SR:i:"):
                SR = int(f[5:])
            if f.startswith("BO:i:"):
                BO = int(f[5:])
                tagged_nodes += 1
            if f.startswith("NO:i:"):
                NO = int(f[5:])
        #Get the length from the sequence if the LN tag is not given
        node[fields[1]] = node[fields[1]]._replace(BO=BO, NO=NO, SO=SO, SR=SR, SN=SN, LN=LN, Seq=fields[2])
        if node[fields[1]].LN == -1:
            node[fields[1]] = node[fields[1]]._replace(LN = len(fields[2]))
    logger.info("\t%d total Nodes Processed"%(total_nodes))
    logger.info("\t%d nodes with tags"%(tagged_nodes))
    
    reader.seek(0)
    while True:
        line = reader.readline()
        if not line:
            break
        if line[0] != 'L':
            continue
        fields = line.split('\t')
        n1 = fields[1]
        o1 = fields[2]
        n2 = fields[3]
        o2 = fields[4]
        if o1 != "+" or o2 != "+" or node[n1].SR != 0 or node[n2].SR != 0:
            continue
        if node[n2].SO != (node[n1].SO + node[n1].LN):
            continue
        new_next = n2
        new_previous = edges[n1].prev
        edges[n1] = edges[n1]._replace(next=new_next, prev=new_previous)
        new_next = edges[n2].next
        new_previous = n1
        edges[n2] = edges[n2]._replace(next=new_next, prev=new_previous)


    reader.close()

def read_assemblies(assembly_list, scaffold_nodes, nodes):
    variants = {}
    bubbles = []
    haplotype_list = []
    contigs = set()
    # Generate the variant INFO (the starting and end node of a bubble)
    for i in range(len(scaffold_nodes) - 1):
        # Take into consideration the chromosome switch
        if nodes[scaffold_nodes[i]].SN != nodes[scaffold_nodes[i+1]].SN:
            continue
        contigs.add(nodes[scaffold_nodes[i]].SN)
        bubble = ">%s>%s"%(scaffold_nodes[i], scaffold_nodes[i+1])
        bubbles.append(bubble)
        variants[bubble] = {}
    logger.info("\nNumber of possible bubbles = %d"%(len(bubbles)))
    gaf_files = []
    with open(assembly_list, 'r') as f:
        for line in f:
            gaf_files.append(line.rstrip())
    for gaf in gaf_files:
        logger.info("\nProcessing Assembly File %s"%(gaf))
        haplotype = gaf.split("/")[-1][:-4]
        haplotype_list.append(haplotype)
        find_variant_alleles(gaf, variants, nodes, haplotype)

    haplotype_list.sort()
    contigs = list(contigs)
    contigs.sort()    
            
    return variants, haplotype_list, contigs

def find_variant_alleles(gaf, variants, nodes, haplotype):
    if is_file_gzipped(gaf):
        reader = gzip.open(gaf, 'rt')
    else:
        reader = open(gaf, 'r')
    counts = defaultdict(lambda: 0)
    while True:
        alignment = reader.readline()
        if not alignment:
            break
        alignment = alignment.split('\t')
        logger.debug("\tProcessing contig %s"%(alignment[0]))
        path = list(filter(None, re.split('(>)|(<)', alignment[5])))
        rv = False
        confused = False
        orient = None
        orient_list = []
        scaffold_index = []
        scaffold_list = []
        counts['Total Number of Alignments'] += 1
        # Reading the entire alignment once to check if the alignment is in the reverse direction.
        # Also storing the scaffold node indices to quickls access the variant bubbles.
        for n, p in enumerate(path):
            if p in ['>', '<']:
                orient = p
                continue
            if nodes[p].NO != 0:
                continue
            scaffold_index.append(n)
            scaffold_list.append(p)
            orient_list.append(orient)
        if len(orient_list) == 0:
            counts['Number of Alignments with no Scaffold Nodes'] += 1
            logger.debug("\t\tSkipping contig since no scaffold nodes were found.")
            continue
        
        # Assign reversed if there are more scaffold nodes with < direction
        # Overall orientation has been reversed
        if orient_list.count('>') < orient_list.count('<'):
            logger.debug('\t\tContig alignment identified as reversed.')
            counts['reversed'] += 1
            scaffold_index.reverse()
            new_scaffold_index = []
            for i in scaffold_index:
                new_scaffold_index.append(len(path) - i)
            scaffold_index = new_scaffold_index
            orient_list.reverse()
            new_orient_list = []
            for n,o in enumerate(orient_list):
                if o == ">":
                    new_orient_list.append("<")
                elif o == "<":
                    new_orient_list.append(">")
            orient_list = new_orient_list
            path = reverse_path(path)
        # Assign confused if there are scaffold nodes in both directions.
        breakpoint = []
        if orient_list.count('>') != 0 and orient_list.count('<') != 0:
            # Break up the path into forward and reverse segments.
            logger.debug('\t\tContig alignment contains inverted scaffold nodes.')
            counts['Number of Alignment with Scaffold Nodes in both orientations'] += 1
            orient = None
            for n,o in enumerate(orient_list):
                if orient == None:
                    orient = o
                    breakpoint.append(orient)
                    continue
                if o != orient:
                    breakpoint.append(n)
                    breakpoint.append(o)
                    orient = o
            breakpoint.append(len(orient_list))
        else:
            breakpoint = [">", len(orient_list)]

        start_path = scaffold_index[0]
        start_index = 0
        for n_bp, b in enumerate(breakpoint):
            if b in ['<', '>']:
                orient  = b
                continue
            end_path = scaffold_index[b-1]
            subpath = path[start_path-1:end_path+1]
            subindex = scaffold_index[start_index:b]
            # If there is only one scaffold node, then I can conisder that to be an inversion event with the inverted allele 
            # Cannot be supported since the inverted node is a scaffold node and hence variant cannot be defined between two nodes that are not consecutive
            # >s1<s2>s3: In this case, since all at scaffold nodes, the variants defined are >s1>s2 and >s2>s3. >s1>s3 is not allowed.
            # Instead now I just ignore these short inversions.
            if len(subindex) == 1:
                logger.debug('\t\tContig alignment contains single scaffold node breakpoint.')
            else:
                diff = subindex[0]
                subindex = [s-diff+1 for s in subindex]
                if orient == '<':
                    subpath = reverse_path(subpath)
                    subindex.reverse()
                    new_scaffold_index = []
                    for i in subindex:
                        new_scaffold_index.append(len(subpath) - i)
                    subindex = new_scaffold_index
                

            s = 0
            e = 1
            while e < len(subindex):
                sn = subpath[subindex[s]]
                en = subpath[subindex[e]]
                bub = ">%s>%s"%(sn,en)
                assert(nodes[sn].NO == 0)
                assert(nodes[en].NO == 0)
                var = subpath[subindex[s]-1:subindex[e]+1]
                try:
                    variants[bub][haplotype].add(''.join(var))
                except KeyError:
                    variants[bub][haplotype] = set()
                    variants[bub][haplotype].add(''.join(var))
                s=e
                e=s+1

            # Updating the breakpoint start variables
            try:
                start_path = scaffold_index[b]
                start_index=b
            except IndexError:
                assert b==len(scaffold_index)
    
    for key, value in variants.items():
        try:
            var = value[haplotype]
            # TODO: Check if the variant bubbles with multiple ATs have the same sequence.
            if len(var) > 1:
                counts['Number of Variant Bubbles with multiple Allele Traversals found'] += 1
                seqs = get_sequences(list(var), nodes)
                seqs = set(seqs)
                if len(seqs) == 1:
                    value[haplotype] = set([list(value[haplotype])[0]])
                    counts['Number of Variant Bubbles with multiple Allele Traversals found but they are the same sequence'] += 1
        except KeyError:
            counts['Number of Variant Bubbles not found in the alignments'] += 1
    for key, value in counts.items():
        logger.info("\t%s: %d"%(key, value))
    reader.close()


def reverse_path(path):
    new_path = []
    for p in path:
        if p == '<':
            orient = ">"
            continue
        elif p == '>':
            orient = "<"
            continue
        new_path.insert(0,p)
        new_path.insert(0,orient)
    return new_path


def write_vcf(writer, variants, ref_alleles, haplotypes, contigs, nodes, keep_all_records, external_variants):
    write_header(writer, contigs, haplotypes)
    write_records(writer, variants, ref_alleles, haplotypes, nodes, keep_all_records, external_variants)
    pass

def write_header(writer, contigs, haplotypes):
    diploid_samples = list(set([x.split('.')[0] for x in haplotypes if '.' in x]))
    diploid_samples.sort()
    haploid_samples = list(set([x.split('.')[0] for x in haplotypes if '.' not in x]))
    haploid_samples.sort()
    samples = diploid_samples+haploid_samples
    
    # Writing basic information
    writer.write("##fileformat=VCFv4.2")
    writer.write("##phased=True")
    
    # Writing INFO Header
    writer.write('##INFO=<ID=CONFLICT,Number=.,Type=String,Description="Assembly names for which there are multiple conflicting allele traversals">')
    writer.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in the assembly panel">')
    writer.write('##INFO=<ID=DIPLOID_AC,Number=A,Type=Integer,Description="Total number of alternate alleles in diploid assemblies">')
    writer.write('##INFO=<ID=HAPLOID_AC,Number=A,Type=Integer,Description="Total number of alternate alleles in haploid assemblies">')
    writer.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1] in the assembly panel">')
    writer.write('##INFO=<ID=DIPLOID_AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1] in diploid assemblies">')
    writer.write('##INFO=<ID=HAPLOID_AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1] in haploid assemblies">')
    writer.write('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data in the assembly panel">')
    writer.write('##INFO=<ID=DIPLOID_NS,Number=1,Type=Integer,Description="Number of samples with data in diploid assemblies">')
    writer.write('##INFO=<ID=HAPLOID_NS,Number=1,Type=Integer,Description="Number of samples with data in haploid assemblies">')
    writer.write('##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">')
    
    # Writing FILTER Header
    writer.write('##FILTER=<ID=PASS,Description="All filters passed">')
    writer.write('##FILTER=<ID=NS80,Description="Number of Diploid Samples with data less than 80%">')
    
    # Writing FOMRAT Header
    writer.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')

    # Writing Contig Header
    for contig in contigs:
        writer.write(contig_header_to_string(contig))
    writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s"%('\t'.join(samples)))


def write_records(writer, variants, ref_alleles, haplotypes, nodes, keep_all_records, external_variants):

    diploid_samples = list(set([x.split('.')[0] for x in haplotypes if '.' in x]))
    diploid_samples.sort()
    haploid_samples = list(set([x.split('.')[0] for x in haplotypes if '.' not in x]))
    haploid_samples.sort()
    num_variants = {'skipped': 0, 'processed': 0}

    ext_var_pointer = 0
    ext_var_chrom = None
    write_begin = False

    for bub in variants.keys():
        
        nd = nodes[bub.split(">")[1]]
        chr = nd.SN
        pos = nd.SO + nd.LN
        if ext_var_chrom != chr or (ext_var_chrom is None and ext_var_pointer == 0):
            ext_var_chrom = chr
            ext_var_pointer = 0
            write_begin = True

        # writing any external variants that are on the starting scaffold of the bubble (only done for the first bubble of the chromosome)
        if write_begin and external_variants and ext_var_pointer < len(external_variants[chr]):
            while (ext_var_pointer < len(external_variants[chr]) and external_variants[chr][ext_var_pointer]['POS'] < pos):
                ext_var = external_variants[chr][ext_var_pointer]
                if ext_var['on_scaffold']:
                    ref = ext_var['REF']
                    alts = ext_var['ALT']
                    id = f'EXT-{chr}-{ext_var["POS"]}'
                    qual = '.'
                    filter = 'PASS'
                    info = '.'
                    format_str = 'GT'
                    genotypes = ['.'] * (len(diploid_samples) + len(haploid_samples))
                    # Writing the external variant record
                    writer.write(f'{chr}\t{ext_var["POS"]}\t{id}\t{ref}\t{",".join(alts)}\t{qual}\t{filter}\t{info}\t{format_str}'+"\t".join(genotypes))
                else:
                    # code block for future use in-case I want to do stuff with the variants in bubbles.
                    pass
                ext_var_pointer += 1
            write_begin = False
        
        variant = variants[bub]
        ref_path = ref_alleles[bub]
        alleles = {ref_path: 0}
        ns = 0
        ns_diploid = 0
        ns_haploid = 0
        count = 1
        allele_count = {0: 0}
        conflict = []
        genotypes = []
              
        # Reading the diploid assembly information
        for sample in diploid_samples:
            #TODO: Hardcoded the identifier for the maternal and paternal haplotypes
            mat = sample+'.2'
            pat = sample+'.1'
            gen = Genotype()
            try:
                mat_allele = variant[mat]
                if len(mat_allele) > 1:
                    mat_allele = ['.']
                    conflict.append(mat)
            except KeyError:
                mat_allele = ['.']
            try:
                pat_allele = variant[pat]
                if len(pat_allele) > 1:
                    pat_allele = ['.']
                    conflict.append(pat)
            except KeyError:
                pat_allele = ['.']
            assert len(mat_allele) == 1 and len(pat_allele) == 1
            mat_allele = list(mat_allele)[0]
            pat_allele = list(pat_allele)[0]
            pat_anum = '.'
            if pat_allele != '.':
                try:
                    pat_anum = alleles[pat_allele]
                    allele_count[pat_anum] += 1
                except KeyError:
                    pat_anum = count
                    alleles[pat_allele] = count
                    count += 1
                    allele_count[pat_anum] = 1
            gen.paternal = pat_anum
            mat_anum = '.'
            if mat_allele != '.':
                try:
                    mat_anum = alleles[mat_allele]
                    allele_count[mat_anum] += 1
                except KeyError:
                    mat_anum = count
                    alleles[mat_allele] = count
                    count += 1
                    allele_count[mat_anum] = 1
            gen.maternal = mat_anum
            genotypes.append(gen.to_string())
            #Determining number of available samples and conflicting samples
            if pat_anum != '.' and mat_anum != '.':
                ns += 1
                ns_diploid += 1
        #Determing filter
        if len(diploid_samples) == 0:
            filter = ['PASS']
        elif ns_diploid/len(diploid_samples) < 0.8:
            filter = ['NS80']
        else:
            filter = ['PASS']
        
        # Reading the haploid haplotype assembly information
        for sample in haploid_samples:
            gen = Haplotype()
            try:
                allele = variant[sample]
                if len(allele) > 1:
                    allele = ['.']
                    conflict.append(sample)
            except KeyError:
                allele = ['.']
            assert len(allele) == 1
            allele = list(allele)[0]
            anum = '.'
            if allele != '.':
                try:
                    anum = alleles[allele]
                    allele_count[anum] += 1
                    gen._remove = True
                except KeyError:
                    anum = count
                    alleles[allele] = count
                    count += 1
                    allele_count[anum] = 1
            gen.allele = anum
            if gen.allele != '.':
                ns += 1
                ns_haploid += 1
            genotypes.append(gen.to_string())
        
        #Determine alternate allele count
        ac = {}
        ac_diploid = {}
        ac_haploid = {}
        for i in range(len(allele_count)):
            ac[i] = 0
            ac_diploid[i] = 0
            ac_haploid[i] = 0
        
        # iterating through diploid assembly genotypes
        for gen in genotypes[0:len(diploid_samples)]:
            g = gen.split('|')
            try:
                ac[int(g[0])] += 1
                ac_diploid[int(g[0])] += 1
            except ValueError:
                pass
            try:
                ac[int(g[1])] += 1
                ac_diploid[int(g[1])] += 1
            except ValueError:
                pass
        
        # iterating through haploidhaplotype genotype
        for gen in genotypes[len(diploid_samples):]:
            try:
                ac[int(gen)] += 1
                ac_haploid[int(gen)] += 1
            except ValueError:
                pass
        
        # determine allele traversals
        at = []
        for key,_ in alleles.items():
            at.append(key)
        seq = get_sequences(at, nodes)
        
        # looking for same alleles with different labels
        allele_to_new = {}
        for i in range(len(seq)):
            allele_to_new[i] = i
        new_seq = []
        for s in seq:
            if s not in new_seq:
                new_seq.append(s)
        # need to make allele traversal list for unique alleles
        at_new = []
        for i in range(len(new_seq)):
            at_new.append(at[seq.index(new_seq[i])])
        
        # remapping the alleles to sequences
        for i in range(len(seq)):
            allele_to_new[i] = new_seq.index(seq[i])
        new_allele_to_sequence = {}
        for old_allele, new_allele in allele_to_new.items():
            if new_allele in new_allele_to_sequence:
                assert seq[old_allele] == new_allele_to_sequence[new_allele]
            else:
                new_allele_to_sequence[new_allele] = seq[old_allele]
        
        # updating the genotypes
        for i, gen in enumerate(genotypes):
            haps = [int(a) if a != '.' else a for a in gen.split('|')]
            new_genotype = []
            for h in haps:
                if h == '.':
                    new_genotype.append('.')
                elif h == 0:
                    new_genotype.append('0')
                else:
                    new_genotype.append(str(allele_to_new[h]))
            genotypes[i] = '|'.join(new_genotype)
        
        # remapping allele counts from old allele to new allele
        new_ac = {}
        new_ac_diploid = {}
        new_ac_haploid = {}
        for old_allele, new_allele in allele_to_new.items():
            try:
                new_ac[new_allele] += ac[old_allele]
                new_ac_diploid[new_allele] += ac_diploid[old_allele]
                new_ac_haploid[new_allele] += ac_haploid[old_allele]
            except KeyError:
                new_ac[new_allele] = ac[old_allele]
                new_ac_diploid[new_allele] = ac_diploid[old_allele]
                new_ac_haploid[new_allele] = ac_haploid[old_allele]

        # preparing allele counts
        new_ac = list(new_ac.values())
        new_ac_diploid = list(new_ac_diploid.values())
        new_ac_haploid = list(new_ac_haploid.values())

        # checking if only one allele (the reference allele) is present.
        # ignore if this is the case
        if len(new_ac) == 1:
            num_variants['skipped'] += 1
            if not keep_all_records:
                continue

        num_variants['processed'] += 1
        
        #Determine allele frequencies
        tot = sum(new_ac)
        af = [x/tot for x in new_ac[1:]]
        tot = sum(new_ac_diploid)
        if tot == 0:
            af_diploid = [0 for _ in new_ac_diploid[1:]]
        else:
            af_diploid = [x/tot for x in new_ac_diploid[1:]]
        tot = sum(new_ac_haploid)
        if tot == 0:
            af_haploid = [0 for _ in new_ac_haploid[1:]]
        else:
            af_haploid = [x/tot for x in new_ac_haploid[1:]]
        
        id = bub
        ref = new_seq[0]
        alt = new_seq[1:]
        # if alt is empty, replace it with a .
        if alt == []:
            alt = ['.']
        qual = 60
        info={"CONFLICT": conflict, "AC": new_ac[1:], "DIPLOID_AC": new_ac_diploid[1:], "HAPLOID_AC": new_ac_haploid[1:], "AF": af, "DIPLOID_AF": af_diploid, "HAPLOID_AF": af_haploid, "NS": ns, "DIPLOID_NS": ns_diploid, "HAPLOID_NS": ns_haploid, "AT": at_new}
        writer.write(variant_record_to_string(chr, pos, id, ref, alt, qual, filter, deepcopy(info), genotypes))
        
        nd = nodes[bub.split(">")[-1]]
        pos = nd.SO + nd.LN
        
        # writing any external variants that are on the starting scaffold of the bubble (only done for the first bubble of the chromosome)
        if external_variants and ext_var_pointer < len(external_variants[chr]):
            while (ext_var_pointer < len(external_variants[chr]) and external_variants[chr][ext_var_pointer]['POS'] < pos):
                ext_var = external_variants[chr][ext_var_pointer]
                if ext_var['on_scaffold']:
                    ref = ext_var['REF']
                    alts = ext_var['ALT']
                    id = f'EXT-{chr}-{ext_var["POS"]}'
                    qual = '.'
                    filter = 'PASS'
                    info = '.'
                    format_str = 'GT'
                    genotypes = ['.'] * (len(diploid_samples) + len(haploid_samples))
                    # Writing the external variant record
                    writer.write(f'{chr}\t{ext_var["POS"]}\t{id}\t{ref}\t{",".join(alts)}\t{qual}\t{filter}\t{info}\t{format_str}\t'+"\t".join(genotypes))
                else:
                    # code block for future use in-case I want to do stuff with the variants in bubbles.
                    pass
                ext_var_pointer += 1

    logger.info("\nNumber of variant records lacking alternate alleles or unavailable alleles: %d", num_variants['skipped'])
    logger.info("Number of variant records in the VCF: %d", num_variants['processed'])
        

def get_reference_alleles(edges, bubbles):
    ref_alleles = {}
    allele = ''
    for bubble in bubbles:
        b = bubble.split(">")
        n1 = b[1]
        n2 = b[-1]
        allele = ''
        allele += '>'+n1
        n = n1
        while n != n2:
            n = edges[n].next
            allele += '>'+n
        ref_alleles[bubble] = allele
    
    return ref_alleles


def get_sequences(at, nodes):
    alleles = []
    for traversal in at:
        path = list(filter(None, re.split('(>)|(<)', traversal)))
        seq = nodes[path[1]].Seq[-1]
        orient = None
        for nd in path[2:-2]:
            if nd in ['>', '<']:
                orient = nd
                continue
            if orient == '>':
                seq += nodes[nd].Seq
            elif orient == '<':
                seq += reverse_complement(nodes[nd].Seq)
        assert seq != ''    
        alleles.append(seq)
    return alleles

def reverse_complement(seq):
    seq = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
    seq = seq.upper()
    # reverse strand
    seq = seq[::-1]
    return seq

def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b'\x1f\x8b'

def write_to_file(line, writer):
    try:
        writer.write(line)
    except TypeError:
        writer.write(str.encode(line))

def format_header_to_string(id, number, type, description):
    return '##FORMAT=<ID=%s,Number=%s,Type=%s,Description="%s">'%(id, number, type, description)

def info_header_to_string(id, number, type, description):
    return '##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">'%(id, number, type, description)

def filter_header_to_string(id, description):
    return '##FILTER=<ID=%s,Description="%s">'%(id, description)

def contig_header_to_string(id):
    return '##contig=<ID=%s>'%(id)

def variant_record_to_string(chr, pos, id, ref, alt, qual, filter, info, genotypes):
    for key, value in info.items():
        if isinstance(value, abc.Iterable):
            info[key] = ','.join([str(a) for a in value])
        else:
            info[key] = str(value)
    return '%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s\tGT\t%s'%(chr, pos-1, id, ref, ','.join(alt), qual, ','.join(filter), ';'.join(['%s=%s'%(key,values) for key, values in info.items()]), '\t'.join(genotypes))

# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    arg('-g', '--gfa', dest='gfa', metavar="GFA", required=True,
        help='GFA file with the sort keys (BO and NO tagged)')
    arg('-dip', '--diploid', dest='diploid', metavar="DIPLOID",
        help='Text file with the list of diploid assembly-to-graph GAF files.')
    arg('-hap', '--haploid', dest='haploid', metavar="HAPLOID",
        help='Text file with the list of haploidhaplotype assembly-to-graph GAF files.')
    arg('-a', '--keep-all-records', dest='keep_all_records', action='store_true',
        help='Flag to keep all the records regardless of whether alt allele is present or all genotypes are unavailable.')
    arg('-e', '--external-vcf', dest='external_vcf', metavar="EXTERNAL_VCF",
        help='External VCF file with small InDels to merge with the generated VCF (The VCF should be sorted).')
    arg('-o', '--output', dest='output', metavar="OUTPUT",
        help='Output VCF path. Default is stdout.')
    
# fmt: on

def validate(args, parser):
    if not args.diploid and not args.haploid:
        parser.error("No assemblies have been provided. Please provide diploid or haploid assemblies.")
    
def main(args):
    run(**vars(args))