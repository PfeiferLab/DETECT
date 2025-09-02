import random
from sys import argv
from Bio.Seq import Seq,MutableSeq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import json
import numpy as np
import math
import argparse
import sys

parser = argparse.ArgumentParser(description="Creates a FASTA file for each haplotype in a VCF.")
parser.add_argument("-i","--input-fasta",dest="input_fasta",help="The input FASTA file.")
parser.add_argument("-v","--input-vcf",dest="input_vcf",help="The input trio VCF file.")
parser.add_argument("-s","--sample-name",dest="sample_name",help="Sample name to turn into FASTA.")
parser.add_argument("-o","--output-stem",dest="output_stem",help="Output stem of FASTAs.")
args = parser.parse_args()

input_file = SeqIO.parse(open(args.input_fasta),'fasta')
input_vcf = open(args.input_vcf)
config = json.load(open('config/config.json'))
var_dict={}
for line in input_vcf:
    if "#" in line and "CHROM" not in line:
        continue

    fields = line.strip().split()
    
    if "CHROM" in line:
        sample_index = fields.index(args.sample_name)
        continue

    chrom = fields[0]
    pos = int(fields[1])
    ref = fields[3]
    alt = fields[4]
    gt1,gt2 = [int(x) for x in fields[sample_index].split('|')]
    if chrom not in var_dict.keys():
        var_dict[chrom] = []
    if ref != "N":    
        var_dict[chrom].append([pos,ref,alt,gt1,gt2])
print("Read in everything")
for item in input_file:
    genome_seq = item.seq
    h1_seq = MutableSeq(genome_seq)
    h2_seq = MutableSeq(genome_seq)
    header = item.id
    for item in var_dict[header]:
        pos = item[0]
        ref = item[1]
        alt = item[2]
        gt1 = item[3]
        gt2 = item[4]
        if gt1 != 0:
            h1_seq[pos-1] = alt
        if gt2 != 0:
            h2_seq[pos-1] = alt

    h1_seq = Seq(h1_seq)
    h2_seq = Seq(h2_seq)

    h1_record = SeqRecord(h1_seq,id=header,description="")
    h2_record = SeqRecord(h2_seq,id=header,description="")
    SeqIO.write(h1_record,args.output_stem+".1.fa","fasta")
    SeqIO.write(h2_record,args.output_stem+".2.fa","fasta")
