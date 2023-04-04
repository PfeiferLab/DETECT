#~/.conda/envs/snakemake/bin/python
import argparse
import sys
import os
import json
import datetime
parser = argparse.ArgumentParser(description="Creates a config file, which is then passed to the DETECT workflow")
required_args = parser.add_argument_group('Required arguments')
required_args.add_argument("-C","--coverages",dest="coverages",help="Comma delimited string with coverages of Dam,Sire,Offspring in empirical data",required = True)
required_args.add_argument("-F","--filter-file",dest="filter_file",help="File that contains filter names, min, max, and stepsize of thresholds.", required = True)
required_args.add_argument("-O","--output-file",dest="output_file",help="Output .txt file.", required = True)
required_args.add_argument("-R","--reference-genome",dest="ref",help="Reference Genome of organism in question", required = True)
required_args.add_argument("-U","--mutation-rate",dest="mutation_rate",help="Estimate mutation rate of organism.", required = True)

optional_args = parser.add_argument_group('Optional arguments')
optional_args.add_argument("-CL","--contig-list",dest="chrom_list",help="List of contigs you want DNMs to be detected on, one per line. Default: All contigs.",default="ALL",required=False)
optional_args.add_argument("-FL","--fragment-length",dest="frag_len",help="Fragment lengths of empirical dataset. Default: 300",default=300,required=False)
optional_args.add_argument("-RL","--read-length", dest="read_length",help="Read lengths of empirical dataset. Default: 100",default=100,required=False)
optional_args.add_argument("-P","--pedigree",dest="pedigree",help="Comma delimited string with names of Dam,Sire,Offspring in VCF. Must be specified if -V is specified. Default: \"parent_1,parent_2,child\"",required=False)
optional_args.add_argument("-SD","--frag-stdv",dest="frag_stdv",help="Fragment length standard deviation. Deafult: 30",default=30,required=False)
optional_args.add_argument("-SP","--software-paths",dest="path_list",help="List of paths to respective bioinformatics softwares. Required softwares are SAMtools, Mason, GATK4, and BWA. Default: Native",required=False)
optional_args.add_argument("-V","--input-variants",dest="input_variants",help="Variant catalog to act as false positives. Optional, but recommended.",default="NONE",required = False)
optional_args.add_argument("--trio",action="store_true",help="Specifies if the input VCF is trio variation data. Either --trio or --population required if -V is specified. Default: None",required = False)
optional_args.add_argument("--population",action="store_true",help="Specifies if the input VCF is population variation data. Either --trio or --population required if -V is specified. Default: None",required = False)
optional_args.add_argument("-WD","--working-directory",dest="working_directory",help="working directory in which to produce files. WARNING, depending on the coverage and size of genome, this workflow may take up significant space. Default: current directory",default=".",required=False)
optional_args.add_argument("--cpus",dest="cpu_count",help="The max number of CPUs you would like to use per job in the workflow. Default: 1",default=1,required=False)
optional_args.add_argument("--num-iterations",dest="num_iterations",help="The number of iterations you would like this to run, to get an idea of the distribution of recommended filters.",default=10,required=False)

args = parser.parse_args()

#Set argparse Variables
cwd = os.path.abspath(os.path.dirname(__file__))
ref = os.path.abspath(args.ref) 
filter_file = open(os.path.abspath(args.filter_file),'r') 
mutation_rate = args.mutation_rate
input_variants = os.path.abspath(args.input_variants)
read_length = args.read_length
read_fragment = args.frag_len
frag_stdv = args.frag_stdv
output_file = os.path.abspath(args.output_file)
coverages = args.coverages 
pedigree = args.pedigree
chrom_list = os.path.abspath(args.chrom_list)
path_list = os.path.abspath(args.path_list)
num_iterations=args.num_iterations
working_directory = os.path.abspath(args.working_directory)

try:
    log_dir = os.path.abspath(args.log_dir)
except:
    log_dir = working_directory
output = {}

#Create Log Directory if not already created
if not os.path.exists(log_dir):
    os.system("mkdir "+log_dir)

#Create Working Directory if not already created
if not os.path.exists(working_directory):
    os.system("mkdir "+working_directory)

#Adding --num-iterations
output['num_iterations'] = num_iterations

# Adding -WD
output["working_directory"] = working_directory

#Adding -R
output["reference_genome"] = ref

#Adding -U
output["mu"] = mutation_rate

#Adding -RL
output["read_length"] = read_length

#Adding -RF
output["read_fragment"] = read_fragment

#Adding -SD
output["frag_stdv"] = frag_stdv

#Adding -O 
output["output_file"] = output_file

if not os.path.exists(input_variants) and  args.input_variants != "NONE":
    print("Input Variant File does not exist!")
    sys.exit()

#Pull from dict file for chrom names/lengths                                    
#Adding "chroms"
dict_file = ref.replace(".fa",".dict")
dict_file = dict_file.replace(".fna",".dict")
chroms = []
if chrom_list != "ALL":
    chrom_list = open(chrom_list,'r')#    for line in chrom_list:
    for line in chrom_list:
        chroms.append(line.strip())

try:
    open_dict_file = open(dict_file,'r')
except:
    print("No dictionary detected. Please rerun this command after creating the dictionary.")                        
    sys.exit()

vcf_idx = input_variants.replace(".vcf",".vcf.idx")

try:
    open_vcf_idx = open(vcf_idx,'r')
except:
    print("No index for input VCF detected. Please rerun this command after indexing you input VCF with gatk IndexFeatureFile.")
    sys.exit()

output["chroms"] = {}
for line in open_dict_file: 
    if len(line.split()) > 2:
        fields = line.strip().split()
        chrom_name = fields[1].split(":")[1]
        chrom_length = fields[2].split(":")[1]
        if chrom_list != "ALL":
            if chrom_name in chroms:
                output["chroms"][chrom_name] = chrom_length
        else:
            output["chroms"][chrom_name] = chrom_length
if len(output["chroms"].keys()) == 0:
    print("No Chromosomes were added. This means your chromosome list does not match any of the chromosomes in your dictionary file.")
    sys.exit()

#Adding -P
pedigree = pedigree.strip().split(",")
output["names"] = {}
output["names"]["parent_1"] = pedigree[0]
output["names"]["parent_2"] = pedigree[1]
output["names"]["child"] = pedigree[2]


#Adding -C                                                                      
coverages = coverages.strip().split(",")                                        
output["coverages"] = {}                                                        
output["coverages"]["parent_1"] = coverages[0]                                 
output["coverages"]["parent_2"] = coverages[1]                                 
output["coverages"]["child"] = coverages[2]

#Adding -F
output["filters"] = {}
for line in filter_file:
    if "#" not in line:
        name,min_filter,max_filter,step=line.strip().split()

        output["filters"][name] = {}
        output["filters"][name]["min"] = str(float(min_filter))
        output["filters"][name]["max"] = str(float(max_filter))
        output["filters"][name]["step"] = str(float(step))

#Adding apps 
output["apps"] = {}
for line in open(path_list):
    fields = line.strip().split()
    output["apps"][fields[0]] = fields[1]

deps = ["bwa","mason_simulator","python","samtools","gatk"]
for dep in deps:
    if dep not in output["apps"].keys():
        output["apps"][dep] = dep

#Adding -V                                                                      
if args.input_variants != "NONE":                                                    
    output["input_variants"] = input_variants                                   

#Setting Trio or Population Level VCF
if args.trio:
    output["trio"] = 1
else:
    output["trio"] = 0

if args.population:
    output["population"] = 1
else:
    output["population"] = 0

if args.trio and args.population:
    print("VCF cannot be both a trio and a population!")
    sys.exit()

if not args.trio and not args.population and args.input_variants != "NONE":
    print("You must specify whether the VCF is a trio or a population!")
    sys.exit()

output["num_cores"] = args.cpu_count

output["snakemake_dir"] = cwd

os.makedirs(working_directory+"/config",exist_ok=True)
outfile = open(working_directory+"/config/config.json",'w')
outfile.write(json.dumps(output,indent=4))
    
print("Config file created in working directory "+working_directory+".")
