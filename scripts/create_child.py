import random
from sys import argv
import json
import numpy as np
import math
import argparse
import sys

parser = argparse.ArgumentParser(description="Creates a VCF file with a synthesized child based on .")
parser.add_argument("-i","--input-vcf",dest="input_vcf",help="The input VCF file.")
parser.add_argument("-o","--output-VCF",dest="output_file",help="Output VCF file.")
args = parser.parse_args()

input_vcf = open(args.input_vcf)
poly_file = open(args.output_file,'w')

ab_mean=0.50
ab_sd=0
total_len = 0
contig_lines=[]

config = json.load(open('config/config.json'))

poly_file.write("##fileformat=VCFv4.2\n")
poly_file.write("##FORMAT=<ID=GT,Number=A,Type=Float,Description=\"Genotype\""+">\n")

chroms = list(config['chroms'].keys())


mut_dict={}

for line in input_vcf:
    if "#" in line and "CHROM" not in line:
        poly_file.write(line)
        continue

    fields = line.strip().split()

    if "CHROM" in line:
        p1_index = fields.index(config['names']['parent_1'])
        p2_index = fields.index(config['names']['parent_2'])
        poly_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+config['names']['parent_1']+"\t"+config['names']['parent_2']+"\t"+config['names']['child']+"\n")
        continue

    chrom = fields[0]
    pos = fields[1]
    ref = fields[3]
    alt = fields[4]
    p1_gt = fields[p1_index]
    p2_gt = fields[p2_index]
    if p1_gt == p2_gt and p1_gt == '0|0':
        continue
    p1_allele = np.random.choice(p1_gt.split('|'),1)[0]
    p2_allele = np.random.choice(p2_gt.split('|'),1)[0]
    ch_gt = p1_allele+'|'+p2_allele
    if chrom in chroms:
        poly_file.write(chrom+'\t'+pos+'\t.\t'+ref+'\t'+alt+'\t'+'1000'+'\t'+'PASS'+'\t'+'.'+'\t'+'GT'+'\t'+p1_gt+'\t'+p2_gt+'\t'+ch_gt+'\n')
