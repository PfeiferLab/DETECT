import argparse

parser = argparse.ArgumentParser(description="Removes contigs from VCF header that are not in the specific chromosomal VCF.")
parser.add_argument("-i","--input-VCF",dest="input_file",help="The input FASTA file.")
parser.add_argument("-c","--chrom",dest="chrom",help = "The chromosome that is present in the VCF.")
parser.add_argument("-o","--output-VCF",dest="output_file",help="Output VCF file.")
args = parser.parse_args()

input_file = open(args.input_file,'r')
output_file = open(args.output_file,'w')
chrom = args.chrom

for line in input_file:
    if '##contig' in line and chrom+"," not in line:
        continue
    else:
        output_file.write(line)
