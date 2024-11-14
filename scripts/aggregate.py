import os
import numpy as np
from sys import argv
import argparse
import json

parser = argparse.ArgumentParser(description="Creates the final output table.")
parser.add_argument("-i","--input-file",dest="input_file",help="The input config file.")
parser.add_argument("-v","--vars-used",dest="vars_used",help="Denoting whether variants were used in the pipeline or not, either a 1 or 0.")
parser.add_argument("-o","--output-file",dest="output_file",help="Output file name.")
parser.add_argument("-m","--multihit-file",dest="multihit_file",help="file that contains the number of multihits in the variants.")
parser.add_argument("-w","--working-directory",dest="working_directory",help="working directory")
parser.add_argument("-n","--num",dest="num",help="num to input into filenames")
args = parser.parse_args()

input_file =json.load(open(args.input_file,'r'))
output_file = open(args.output_file,'w')
complete_file = open("complete_file.txt",'w')
vars_used = int(args.vars_used)
num=str(args.num)
num_multihits = int(open(args.multihit_file).readlines()[0])
output = []
wd = args.working_directory
output_file.write("Filter\tThreshold\tTotal_Sites\tMutations\tVariants\tSequencing_Errors\n")
complete_file.write("Filter\tThreshold\tTotal_Sites\tMutations\tVariants\tSequencing_Errors\n")

if vars_used:
    orig_var_sites = int(os.popen("grep -v '#' "+wd+"/pipeline/mutations/polymorphisms."+num+".vcf | wc -l").read())
    mv_var_sites = int(os.popen("grep -v '#' "+wd+"/pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV."+num+".polymorphisms.vcf | wc -l").read())
else:
    orig_var_sites = 0
    mv_var_sites = 0
orig_mut_sites = int(os.popen("grep -v '#' "+wd+"/pipeline/mutations/mutations."+num+".vcf | wc -l").read())
mv_mut_sites = int(os.popen("grep -v '#' "+wd+"/pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV."+num+".mutations.vcf | wc -l").read())

original_sites = orig_mut_sites + orig_var_sites
mv_sites = int(os.popen("grep -v '#' "+wd+"/pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV."+num+".vcf | wc -l").read())
mv_errs = mv_sites - mv_var_sites - mv_mut_sites 
output_file.write("original"+'\t'+"N/A"+'\t'+str(original_sites)+'\t'+str(orig_mut_sites)+'\t'+str(orig_var_sites)+'\t'+"0"+'\n')
output_file.write("MVs"+'\t'+"N/A"+'\t'+str(mv_sites)+'\t'+str(mv_mut_sites)+'\t'+str(mv_var_sites)+'\t'+str(mv_errs)+'\n')
complete_file.write("original"+'\t'+"N/A"+'\t'+str(original_sites)+'\t'+str(orig_mut_sites)+'\t'+str(orig_var_sites)+'\t'+"0"+'\n')
complete_file.write("MVs"+'\t'+"N/A"+'\t'+str(mv_sites)+'\t'+str(mv_mut_sites)+'\t'+str(mv_var_sites)+'\t'+str(mv_errs)+'\n')
best_filters = {}


for filter_name in list(input_file['filters'].keys()):
    max_value=input_file['filters'][filter_name]['max']
    min_value=input_file['filters'][filter_name]['min']
    step=input_file['filters'][filter_name]['step']
    for val in np.round(np.arange(float(min_value),float(max_value) + float(step) / 2,float(step)),2):
        total_sites = int(os.popen("grep -v '#' "+wd+"/pipeline/filters/"+filter_name+"/all_chr_trio.downsampled.sorted.mark_dups.MV."+filter_name+"."+str(val)+"."+num+".vcf | wc -l").read())
        mut_sites = int(os.popen("grep -v '#' "+wd+"/pipeline/filters/"+filter_name+"/all_chr_trio.downsampled.sorted.mark_dups.MV."+filter_name+"."+str(val)+"."+num+".mutations.vcf | wc -l").read())
        
        if vars_used:
            var_sites = int(os.popen("grep -v '#' "+wd+"/pipeline/filters/"+filter_name+"/all_chr_trio.downsampled.sorted.mark_dups.MV."+filter_name+"."+str(val)+"."+num+".polymorphisms.vcf | wc -l").read())
        
        else:
            var_sites = 0
        errors = total_sites - mut_sites - var_sites
        complete_file.write(filter_name+"\t"+str(val)+"\t"+str(total_sites).strip()+"\t"+str(mut_sites).strip()+"\t"+str(var_sites).strip()+"\t"+str(errors)+'\n')#debug line
        
        if filter_name not in best_filters.keys():
            best_filters[filter_name] = [val,total_sites,mut_sites,var_sites,errors]
        
        if mut_sites > best_filters[filter_name][2]:
            best_filters[filter_name] = [val,total_sites,mut_sites,var_sites,errors]
        if mut_sites >= best_filters[filter_name][2] and var_sites < best_filters[filter_name][3]:
            best_filters[filter_name] = [val,total_sites,mut_sites,var_sites,errors]
best_filters_list = list(best_filters.keys())

for item in best_filters_list:
    output_file.write(item+'\t'+str(best_filters[item][0])+'\t'+str(best_filters[item][1])+'\t'+str(best_filters[item][2])+'\t'+str(best_filters[item][3])+'\t'+str(best_filters[item][4])+'\n')
    complete_file.write(item+'\t'+str(best_filters[item][0])+'\t'+str(best_filters[item][1])+'\t'+str(best_filters[item][2])+'\t'+str(best_filters[item][3])+'\t'+str(best_filters[item][4])+'\n')
output_file.write("Num Multihit Sites:"+str(num_multihits))
complete_file.write("Num Multihit Sites:"+str(num_multihits))
