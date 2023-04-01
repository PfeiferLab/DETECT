import random
from sys import argv
from Bio import Seq
from Bio import SeqIO
import numpy as np
import math
import scipy.stats as st
import argparse
import sys
seq = ['A','G','C','T']
num_list = [0,1,2,3]

###############################################################
####Mean and SD are from Tatsumoto nonMV, filtered variants####
####mean=0.4945519#############################################
####sd=0.05927348##############################################
###############################################################


parser = argparse.ArgumentParser(description="Creates a VCF file with a mutation rate equivalent to the mu input in positions bounded by the input FASTA.")
parser.add_argument("-i","--input-fasta",dest="input_file",help="The input FASTA file.")
parser.add_argument("-v","--input-vcf",dest="input_vcf",help="",default=[],required=False)
parser.add_argument("-u","--mutation-rate",dest="mutation_rate",help="The input FASTA file.")
parser.add_argument("-o","--output-VCF",dest="output_file",help="Output VCF file.")
parser.add_argument("-p","--pedigree",dest="pedigree",help="comma delimited string for dam,sire,offspring")
parser.add_argument("-c","--chromosomes",dest="chroms",help="comma delimited string of chromosomes desired to be mutated.")
parser.add_argument("-n","--num",dest="num",help="the number to add to the end.")
args = parser.parse_args()

input_file = SeqIO.parse(open(args.input_file),'fasta')
mutation_rate = float(args.mutation_rate)
mutation_file = open(args.output_file,'w')
pedigree_list = args.pedigree.strip().split(",")
chrom_list = args.chroms.strip().split(",")
input_vcf = args.input_vcf
num=args.num

ab_mean=0.50
ab_sd=0
total_len = 0
contig_lines=[]
mutation_lines=[]

present_vars={}
if type(input_vcf) == str and input_vcf != "NONE":
    #print("Entered")
    vcf_file = open(input_vcf,'r') 
    for line in vcf_file:
        if "CHROM" in line:
            #print("Entered")
            names=line.strip().split()[-3:]
            for i in range(0,len(names)):
                if names[i] == pedigree_list[0]:
                    p1_index = i
                elif names[i] == pedigree_list[1]:
                    p2_index = i
                elif names[i] == pedigree_list[2]:
                    ch_index = i
                else:
                    print("ERROR. Names in pedigree list do not match names in VCF!")
        if "#" not in line:
            fields=line.split()
            chrom=fields[0]
            pos=int(fields[1])
            ref_allele=fields[3]
            alt_allele=fields[4]
            p1_gt=fields[p1_index-3].split(":")[0]#.replace("|","/")
            p2_gt=fields[p2_index-3].split(":")[0]#.replace("|","/")
            ch_gt=fields[ch_index-3].split(":")[0]#.replace("|","/")
            present_vars[(chrom,pos)]=[[ref_allele,alt_allele],[p1_gt,p2_gt,ch_gt],0]
#print(present_vars[("test_1","74")])
#print("p1_index:"+str(p1_index))
#print("p2_index:"+str(p2_index))
#print("ch_index:"+str(ch_index))

multihit_count=0
for item in input_file:
    genome_1_seq = item.seq
    header = item.id
    if header in chrom_list:
        mut_pos = []
        total_len+=len(genome_1_seq)
        mean = 2 * mutation_rate * len(genome_1_seq)
        print("mean",mean)
        num_mutations = np.random.poisson(lam=mean)
       
        #if random.random() < 0.5:
        #    num_mutations = int(math.ceil(random.normalvariate(mean,math.sqrt(mean))))
        #else:
        #    num_mutations = int(math.floor(random.normalvariate(mean,math.sqrt(mean))))
        print("num_mutations:"+str(num_mutations))
        mutation = ""
        contig_lines.append("##contig=<ID="+header+",length="+str(len(genome_1_seq))+'>\n')
        mutations = random.sample(range(1,len(genome_1_seq)+1),k=num_mutations)
        mutations.sort()
    #    print(mutations)
        for pos in mutations:
            mutation = ""
            #input((header,pos))
            if (header,pos) in present_vars.keys():
                multihit_count+=1
                print("MULTIHIT:"+str(header)+" "+str(pos))#str(present_vars[(header,pos)]))
                ref = present_vars[(header,pos)][0][0]
                alt = present_vars[(header,pos)][0][1]
                ch_gt = present_vars[(header,pos)][1][2]
                mutation = ""
                if ch_gt == "0/0" or ch_gt == "0|0":
                #    print("ch is homref")
                    while mutation == ref or mutation == "":
                        mutation = random.choice(seq)
                    if mutation == alt:
                        present_vars[(header,pos)][1][2] = random.choice(["0|1","1|0"])
                    else:
                        present_vars[(header,pos)][0][1]+=","+mutation
                        present_vars[(header,pos)][1][2] = random.choice(["0|2","2|0"])
                elif ch_gt == "0/1" or ch_gt == "0|1" or ch_gt == "1|0":
               #     print("ch is het")
                    mutation = random.choice(seq)
                    if mutation == alt:
                        present_vars[(header,pos)][1][2] = "1|1"
                    elif mutation == ref:
                        present_vars[(header,pos)][1][2] = "0|0"
                    else:
                        mutant_gt = random.choice(["0|2","2|0","1|2","2|1"])
                #        print(mutant_gt+" was chosen")
                        present_vars[(header,pos)][0][1]+=","+mutation
                        present_vars[(header,pos)][1][2] = mutant_gt
                elif ch_gt == "1/1" or ch_gt == "1|1":
              #      print("ch is homalt")
                    while mutation == alt or mutation == "":
                        mutation = random.choice(seq)
                    if mutation == ref:
                        present_vars[(header,pos)][1][2] = random.choice(["0|1","1|0"])
                    else:
                        present_vars[(header,pos)][0][1]+=","+mutation
                        present_vars[(header,pos)][1][2] = random.choice(["1|2","2|1"])
                else:
                    print("PROBLEM. CH_GT IS NOT BIALLELIC GENOTYPE YOU WOULD EXPECT")
                    sys.exit()
                #final_gt = '\t'.join(present_vars[(header,pos)][1])
                #final_gt=present_vars[(header,pos)][1][0]+'\t'+present_vars[(header,pos)][1][1]+'\t'+present_vars[(header,pos)][1][2]
             #   print("mutation:"+mutation)
                #reference_allele = ref
                #alt = present_vars[(header,pos)][0][1]
            #    print("final line:"+str(header)+"\t"+str(pos+1)+"\t.\t"+reference_allele+"\t"+alt+"\t1000\tPASS\tMT=1\tGT\t"+final_gt)
                present_vars[(header,pos)][2] = 1
            else:
                while mutation == "" or mutation.upper() == genome_1_seq[pos-1].upper():
                    mutation = random.choice(seq)
                alt = mutation
                reference_allele = genome_1_seq[pos-1].upper()
                present_vars[(header,pos)] = [[reference_allele,alt],["0|0","0|0",random.choice(["0|1","1|0"])],1]
                #present_vars[(header,pos)][1][0] = "0/0"
                #present_vars[(header,pos)][1][1] = "0/0"
                #present_vars[(header,pos)][1][2] = "0/1"
            #present_vars[(header,pos)][2] = 1
            #mutation_lines.append(str(header)+"\t"+str(pos+1)+"\t.\t"+reference_allele+"\t"+alt+"\t1000\tPASS\tAF="+str(final_gt.count("1")/6.0)+";MT=1\tGT\t"+final_gt+"\n")
            #mutation_lines.append(str(header)+"\t"+str(pos+1)+"\t.\t"+reference_allele+"\t"+alt+"\t1000\tPASS\tMT=1\tGT\t"+final_gt+"\n")

#print(str(genome_1_seq))

mutation_file.write("##fileformat=VCFv4.2\n")                               
#mutation_file.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency among genotypes, for each ALT allele, in the same order as listed\">\n")
mutation_file.write("##INFO=<ID=ST,Number=A,Type=Float,Description=\"input mu:"+str(mutation_rate)+"\tmu_calculation:"+str(float(num_mutations)/total_len)+"\">\n")
mutation_file.write("##INFO=<ID=MT,Number=A,Type=Integer,Description=\"Whether the site is a mutation or not; 1: Yes, 0: No\">\n")
mutation_file.write("##FORMAT=<ID=GT,Number=A,Type=Float,Description=\"Genotype\""+">\n")
for line in contig_lines:
    mutation_file.write(line)
mutation_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+pedigree_list[0]+"\t"+pedigree_list[1]+"\t"+pedigree_list[2]+"\n")

for key in sorted(present_vars.keys(),key=lambda element: (element[0],element[1])): 
    ref = present_vars[key][0][0]
    alt = present_vars[key][0][1]
    mut = present_vars[key][2]
    gt = '\t'.join(present_vars[key][1])
    mutation_file.write(key[0]+'\t'+str(key[1])+'\t'+'.'+'\t'+ref+'\t'+alt+'\t'+'1000'+'\t'+'PASS'+'\t'+'MT='+str(mut)+'\t'+'GT'+'\t'+gt+'\n')
print("mulihits:"+str(multihit_count))

                                                            
multihit_file = open('pipeline/mutations/multihit_count.'+str(num)+'.txt','w')
multihit_file.write(str(multihit_count)) 
                                                                         


