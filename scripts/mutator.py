import random
from sys import argv
from Bio import Seq
from Bio import SeqIO
import json
import numpy as np
import math
import argparse
import sys
seq = ['A','G','C','T']
num_list = [0,1,2,3]

parser = argparse.ArgumentParser(description="Creates a VCF file with a mutation rate equivalent to the mu input in positions bounded by the input FASTA.")
parser.add_argument("-i","--input-fasta",dest="input_file",help="The input FASTA file.")
parser.add_argument("-u","--mutation-input",dest="mutation_input",help="The mutation rate/count/file.")
parser.add_argument("-o","--output-VCF",dest="output_file",help="Output VCF file.")
parser.add_argument("-p","--pedigree",dest="pedigree",help="comma delimited string for dam,sire,offspring")
parser.add_argument("-n","--num",dest="num",help="the number to add to the end.")
args = parser.parse_args()

input_file = SeqIO.parse(open(args.input_file),'fasta')
mutation_input = args.mutation_input
mutation_file = open(args.output_file,'w')
pedigree_list = args.pedigree.strip().split(",")
num=args.num

ab_mean=0.50
ab_sd=0
total_len = 0
contig_lines=[]

config = json.load(open('config/config.json'))

mutation_file.write("##fileformat=VCFv4.2\n")
#mutation_file.write("##INFO=<ID=ST,Number=A,Type=Float,Description=\"input mu:"+str()+"\tmu_calculation:"+str(float(num_mutations)/total_length)+"\">\n")
mutation_file.write("##INFO=<ID=MT,Number=A,Type=Integer,Description=\"Whether the site is a mutation or not; 1: Yes, 0: No\">\n")
mutation_file.write("##FORMAT=<ID=GT,Number=A,Type=Float,Description=\"Genotype\""+">\n")

mutation_input_type = ""

chroms = list(config['chroms'].keys())

mutation_input_type = config['mutation_input_type']

mut_dict={}

if mutation_input_type!="file":
    verb="created"
    chrom_lengths = [int(x) for x in config['chroms'].values()]
    total_length = sum(chrom_lengths)
    weights = [x/total_length for x in chrom_lengths]
    
    if mutation_input_type == "count":
        num_mutations = int(mutation_input)
    
    if mutation_input_type == "rate":
        expected_number = float(mutation_input) * total_length * 2
        num_mutations = np.random.poisson(expected_number)
    
    mut_chrom_list = np.random.choice(chroms,size=int(num_mutations),replace=True, p=weights)
    for chrom in mut_chrom_list:
        all_pos = []
        chrom_length = int(config['chroms'][chrom])
        if chrom not in mut_dict.keys():
            mut_dict[chrom] = []
        pos=-1 #initiate pos before while loop, but make it a position that can't normally exist
        while pos in mut_dict[chrom] or pos < 0:
            pos = np.random.randint(1,chrom_length+1)
        mut_dict[chrom].append(pos)

else:
    verb="imported"
    num_mutations=0
    mutation_input_file = open(mutation_input)
    for line in mutation_input_file:
        if "#" in line:
            continue
        fields = line.strip().split()
        chrom=fields[0]
        pos=int(fields[1])
        if chrom not in mut_dict.keys():
            mut_dict[chrom] = []
        mut_dict[chrom].append(pos)
        num_mutations+=1

mutation_lines=[]
print(str(num_mutations)+" mutations "+verb+".")
for item in input_file:
    genome_seq = item.seq
    header = item.id
    if header in mut_dict.keys():
        mut_dict[header].sort()
        mutation_file.write("##contig=<ID="+header+",length="+str(len(genome_seq))+'>\n')
        for pos in mut_dict[header]:
            mutation = ""
            ref = genome_seq[pos-1].upper()
            while mutation == ref or mutation == "":
                mutation = random.choice(seq)
            mutation_lines.append(header+'\t'+str(pos)+'\t.\t'+ref+'\t'+mutation+'\t'+'1000'+'\t'+'PASS'+'\t'+'MT=1'+'\t'+'GT'+'\t'+'0|0'+'\t'+'0|0'+'\t'+random.choice(["0|1","1|0"])+'\n')

mutation_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+pedigree_list[0]+"\t"+pedigree_list[1]+"\t"+pedigree_list[2]+"\n")

for item in mutation_lines:
    mutation_file.write(item)
