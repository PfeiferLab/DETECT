#!/usr/bin/env snakemake
import os
import numpy as np
cwd = os.path.abspath(os.path.dirname(__file__))
SNAKEDIR = config["snakemake_dir"]
workdir: config["working_directory"]
tmpdir: config["working_directory"]
num_chroms = len(config["chroms"].keys())
read_cutoff = 1e6
#TODO Clean up shell commands


import pprint
import numpy
from numpy import random
pp = pprint.PrettyPrinter(indent=4)

wildcard_constraints:
    indiv="|".join(list(config['names'].keys())).replace("_", "\_"), 
    chromosome="|".join(list(config['chroms'].keys())).replace("_", "\_"),
    genome="0|1",
    iter='|'.join([str(x) for x in numpy.arange(1,9999)]),
    run="|".join([str(x) for x in numpy.arange(1,9999)])

#def get_all_input():
#    target_list = ['pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.{iter}.mutations.vcf']
#    filters = get_filters()
#    for filter in filters:
#    	for iter in range(1,int(config['num_iterations'])+1):
#        	target_list += expand('pipeline/filters/{filter}/all_chr_trio.downsampled.sorted.mark_dups.MV.{filter}.{filter_threshold}.{{iter}}.mutations.vcf', filter_threshold = numpy.round(numpy.arange(float(config["filters"][filter]["min"]),float(config["filters"][filter]["max"])+float(config["filters"][filter]["step"])/2,float(config["filters"][filter]["step"])),2),filter=filter),
#        	if 'input_variants' in config.keys():
#            		target_list += expand('pipeline/filters/{filter}/all_chr_trio.downsampled.sorted.mark_dups.MV.{filter}.{filter_threshold}.{{iter}}.polymorphisms.vcf', filter_threshold = numpy.round(numpy.arange(float(config["filters"][filter]["min"]),float(config["filters"][filter]["max"])+float(config["filters"][filter]["step"])/2,float(config["filters"][filter]["step"])),2),filter=filter)
#    if 'input_variants' in config.keys():
#        target_list.append('pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.{iter}.polymorphisms.vcf')
#    return target_list

#def get_all_super_input():
#	target_list = []
#	target_list += expand('pipeline/output/final_{iter}',iter = range(1,int(config['num_iterations'])+1))
#	return ",".join(target_list)

def get_Mutate_input(wildcards):
	if config["trio"]:
		return f'pipeline/input_variants/input_trio.renamed.vcf'
	elif config['population']:
		return f'pipeline/simulated_vcf/simulated_trio_variants.vcf'

def format_vcfmerge(wildcards,input):
	if input.vcf != []:
		return f'-v '+input.vcf
	else:
		return ''

def format_simulate_variants(wildcards,input):
	if input.vcf != []:
		return f'-iv '+input.vcf
	else:
		return ''

def get_mem_mb(wildcards,attempt,input):
	return 64000 #min(max(attempt*(150 *1024),input.size_mb* 1024),307200)

def format_input(wildcards, input):
    return [f'-I {a}' for a in input]

def format_merge_input(wildcards, input):
	return[f'{a}' for a in input]

def format_samtofastq_input(wildcards):
	if wildcards.indiv == 'child':
		return 'pipeline/reads/child_{iter}.golden.mutated.bam'
	else:
		return 'pipeline/reads/{indiv}_{iter}.golden.bam'

def format_genotype(wildcards,input):
    return [f'--variant {a}' for a in input]

def get_read_count(wildcards): 
    read_length = float(config["read_length"])
    chrom_length = float(sum([int(x) for x in config["chroms"].values()]))
    coverage = float(config["coverages"][wildcards.indiv])
    return int(round(((chrom_length*coverage)/read_length)*1.2)/2) 


def get_reads_per_run(wildcards):
    read_count = get_read_count(wildcards)
    num_runs = 16
    num_reads = int(read_count/num_runs)
    return num_reads

rule all:
    input:
        config['output_file']

rule RenameInputTrioVCF:
	input:
		input_file = lambda x: config['input_variants'] if 'input_variants' in config.keys() else 0
	output:
		'pipeline/input_variants/input_trio.renamed.vcf'
	params:
		p1 = config["names"]["parent_1"],
		p2 = config["names"]["parent_2"],
		ch = config["names"]["child"],
		python = config["apps"]["python"],
		snakedir = SNAKEDIR
	shell:
		'{params.python} {params.snakedir}/scripts/rename_samples.py {input.input_file}  \"{params.p1},{params.p2},{params.ch}\" > {output} && gatk IndexFeatureFile -I {output}'

rule Mutate:                                                                    
    input:                                                                      
        fa = config['reference_genome'],
		
    output:                                                                     
        muts='pipeline/mutations/input_mutations.{iter}.phased.vcf',
#	multihit='pipeline/mutations/multihit_count.{iter}.txt'
    params:                                                                     
    	mu  = config['dnm_file'] if config['dnm_file'] != "False" else config['mu'],
	script = 'mutator.py' if config['dnm_file'] != 'False'	else 'mutator_no_dnm_list.py',
        python = config["apps"]["python"],                                      
        chroms = ",".join(list(config["chroms"].keys())),                       
        #vcfmerge = format_vcfmerge, 
        snakedir = SNAKEDIR,
	gatk = config["apps"]["gatk"]                                                     
    resources:
        mem_mb = 8000,
	runtime='15m'
    shell:                                                                      
        '{params.python} {params.snakedir}/scripts/{params.script} -i {input.fa} -u {params.mu} -p \"parent_1,parent_2,child\" -c \"{params.chroms}\" -n {wildcards.iter} -o {output.muts}; '
	'{params.gatk} IndexFeatureFile -I {output.muts}'
#rule PhaseMutations:
#	input:
#		'pipeline/mutations/all_mutations.{iter}.vcf'
#	output:
#		'pipeline/mutations/all_mutations.{iter}.phased.vcf'
#	params:
#		snakedir=SNAKEDIR
#	conda: 'DETECT_053124_env_2.yml'
#	shell:
#		'beagle ped={params.snakedir}/trio.ped gt={input} out=pipeline/mutations/all_mutations.{iter}.phased && gunzip pipeline/mutations/all_mutations.{iter}.phased.vcf.gz'

#rule FilterMutations:
#    input:
#        'pipeline/mutations/all_mutations.{iter}.phased.vcf'
#    output:
#        'pipeline/mutations/mutations.{iter}.phased.vcf'
#    params:
#        gatk = config["apps"]["gatk"]
#    resources:
#        runtime='15m'
#    retries: 20
#    shell:
#       '{params.gatk} SelectVariants -V {input} --select "MT==1" -O {output}'

rule ReformatMutations:
	input:
		'pipeline/mutations/input_mutations.{iter}.phased.vcf'
	output:
		'pipeline/mutations/input_mutations.{iter}.phased.reformatted.vcf'
	params:
		snakedir=config["snakemake_dir"],
		python=config['apps']['python'],
		gatk=config['apps']['gatk']
	resources:
		runtime='15m'
	shell:
		'{params.python} {params.snakedir}/scripts/reformat_vcf.py {input} > {output} && gatk IndexFeatureFile -I {input}'

#rule FilterPolymorphisms:
#    input:
#        'pipeline/mutations/all_mutations.{iter}.phased.vcf'
#    params:
#        gatk = config["apps"]["gatk"]
#    output:
#        'pipeline/mutations/polymorphisms.{iter}.phased.vcf'
#    resources:
#        runtime='15m'
#    shell:
#        '{params.gatk} SelectVariants -V {input} --select "MT==0" -O {output}'

#rule SplitInputMutations:
#	input:
#		input_file = 'pipeline/mutations/all_mutations.{iter}.phased.vcf',
#	output:
#		'pipeline/mutations/{chromosome}.mutations.{iter}.vcf'
#	params:
#		gatk = config["apps"]["gatk"]
#	resources:
#		runtime='15m'
#	shell:
#		'{params.gatk} SelectVariants -V {input.input_file} -L {wildcards.chromosome} -O {output}'

rule splitVCF:
	input:
		muts='pipeline/input_variants/input_trio.renamed.vcf'#config['input_variants'] #pipeline/mutations/all_mutations.{iter}.phased.vcf'
	output:
		'pipeline/mutations/polymorphisms.{chromosome}.phased.vcf'
	resources:
		runtime='15m'
	retries: 20
	shell:
		'gatk SelectVariants -V {input} -L {wildcards.chromosome} -O {output} && gatk IndexFeatureFile -I {output}'

rule VCFtoFasta:
	input:
		ref=config['reference_genome'],
		vcf='pipeline/mutations/polymorphisms.{chromosome}.phased.vcf'
	output:
	#2_parent_1_chr8:1.fa 
		'pipeline/ref/{indiv}_{chromosome}:{genome}.fasta' #expand('pipeline/ref/{{iter}}_{indiv}_{{chromosome}}:{genome}.fa',indiv=config['names'],genome=[0,1]) #format_split_vcf_fasta 'pipeline/ref/{indiv}_{iter}_{genome}.fa'#Akira_tm_autosomal.snp.noMV.recode.pass.0.fa'
	params:
		ref = config['reference_genome'],
		
	#conda: 'DETECT_053124_env.yml'
	resources:
		runtime='15m'
	shell:
		'mkdir -p pipeline/ref && cd pipeline/ref && vcf2fasta -f {input.ref} ../mutations/polymorphisms.{wildcards.chromosome}.phased.vcf'

rule MergeFastas:
	input:
		expand('pipeline/ref/{{indiv}}_{chromosome}:{{genome}}.fasta',chromosome=config['chroms'].keys())
	output:
		'pipeline/ref/{indiv}_renamed.{genome}.fa'
	shell:
		'cat {input} | sed "s/:{wildcards.genome}//g; s/{wildcards.indiv}_//g" > {output}'

rule SimulateReads:
    input:
    	fa = 'pipeline/ref/{indiv}_renamed.{genome}.fa',
    output:
        r1 = temp('pipeline/reads/{indiv}_{run}_{iter}_{genome}.golden.R1.fq'), 
        r2 = temp('pipeline/reads/{indiv}_{run}_{iter}_{genome}.golden.R2.fq'),
	golden_bam = 'pipeline/reads/{indiv}_{run}_{iter}_{genome}.golden.bam'
    params:
        read_length = config['read_length'],
        read_fragmean = config['read_fragment'],
        mason = config["apps"]["mason_simulator"],
        read_count = get_reads_per_run,
    retries: 20
    resources:
    	mem_mb=get_mem_mb,
	runtime='4h'
    threads:
        int(config['num_cores'])
    shell:
        '{params.mason} --read-name-prefix "simulated_{wildcards.indiv}_{wildcards.run}_{wildcards.iter}_{wildcards.genome}" --seed $RANDOM --num-threads {threads} -ir {input.fa} --fragment-mean-size {params.read_fragmean} --illumina-read-length {params.read_length} -n {params.read_count} -o {output.r1} -or {output.r2} -oa {output.golden_bam}' 

rule MergeGoldenBam:
	input:
		expand('pipeline/reads/{{indiv}}_{run}_{{iter}}_{genome}.golden.bam',run=numpy.arange(1,9,1),genome=[0,1])
	params:
		input_list = format_merge_input,
		samtools=config['apps']['samtools']
		#gatk=config['apps']['gatk']
	output:
		'pipeline/reads/{indiv}_{iter}.golden.bam'
	shell:
		'{params.samtools} merge {output} {params.input_list}'
		#'{params.gatk} MergeSamFiles {params.input_list} -O {output}'

rule MutateGoldenBam:
	input:
		bam='pipeline/reads/child_{iter}.golden.bam',
		vcf='pipeline/mutations/input_mutations.{iter}.phased.reformatted.vcf'
	output:
		'pipeline/reads/child_{iter}.golden.mutated.bam'
	params:
		gatk=config['apps']['gatk']
	shell:
		'java -jar ~/jvarkit/dist/biostar404363.jar -p {input.vcf} -o {output} {input.bam}'

rule SamToFastq:
	input:
		format_samtofastq_input
	output:
		fq1='pipeline/reads/{indiv}_{iter}.R1.fq',
		fq2='pipeline/reads/{indiv}_{iter}.R2.fq'
	params:
		gatk=config['apps']['gatk']
	shell:
		'{params.gatk} SamToFastq -I {input} -F {output.fq1} -F2 {output.fq2}'

rule MapReads: 
	input:
		reference = config['reference_genome'],
		fq1 = 'pipeline/reads/{indiv}_{iter}.R1.fq', 
		fq2 = 'pipeline/reads/{indiv}_{iter}.R2.fq',
	threads:
		int(config['num_cores'])
	output:
		'pipeline/sorted_bams/{indiv}.{iter}.bam'
	params:
		bwa = config["apps"]["bwa"],
		samtools = config['apps']['samtools']
	resources:
		mem_mb = 2000 * int(config['num_cores']),
		runtime='12h'
	shell:
		'{params.bwa} mem -t {threads} -R \"@RG\\tID:{wildcards.indiv}_{wildcards.iter}\\tSM:{wildcards.indiv}\\tPL:ILLUMINA\" {input.reference}  {input.fq1} {input.fq2} | samtools view -b - > {output} '

rule DownsampleBam:
	input:
		'pipeline/sorted_bams/{indiv}.{iter}.bam'
	output:
		'pipeline/sorted_bams/{indiv}.downsampled.{iter}.bam'
	params:
		samtools = config['apps']['samtools']
	threads:
		int(config['num_cores'])
	resources:
		mem_mb = 2000 * int(config['num_cores']),
		runtime='4h'
	
	shell:
		'{params.samtools} view -@ {threads} -b -s 0.83333333 {input} > {output}'

rule SortBam:
    input:
        'pipeline/sorted_bams/{indiv}.downsampled.{iter}.bam'
    output:
        'pipeline/sorted_bams/{indiv}.downsampled.sorted.{iter}.bam'
    params:
        gatk = config['apps']['gatk'],
    resources:
        mem_mb=get_mem_mb
    threads: 12
    shell:
        'gatk SortSam --java-options "-Xmx{resources.mem_mb}m" -I {input} --TMP_DIR pipeline/sorted_bams/ -O {output} -SO coordinate '
	

rule MarkDuplicates: 
    input:
        'pipeline/sorted_bams/{indiv}.downsampled.sorted.{iter}.bam'
    output:
        output_bam = 'pipeline/mark_dups/{indiv}.downsampled.sorted.mark_dups.{iter}.bam',
        output_metrics = 'pipeline/mark_dups/{indiv}.downsampled.sorted.mark_dups.metrics.{iter}.txt'
    params:
        gatk = config['apps']['gatk']
    retries: 20
    resources:
    	runtime='1d',
        mem_mb=get_mem_mb,
    shell:
        '{params.gatk} --java-options "-Xmx{resources.mem_mb}m" MarkDuplicates -I {input} -O {output.output_bam} -M {output.output_metrics}  --TMP_DIR pipeline/mark_dups/ --CREATE_INDEX true'

rule BQSR:
	input:
		bam = lambda x: 'pipeline/mark_dups/{indiv}.downsampled.sorted.mark_dups.{iter}.bam' if 'input_variants' in config.keys() else [],
		known_variants = lambda x: config['input_variants'] if 'input_variants' in config.keys() else [],
		reference = config['reference_genome']
	output:
		'pipeline/BQSR/{indiv}.downsampled.sorted.mark_dups.BQSR.{iter}.txt'
	params:
		gatk = config['apps']['gatk']
	resources:
		tmpdir='pipeline/BQSR/'
	shell:
		'{params.gatk} BaseRecalibrator -R {input.reference} -I {input.bam} --known-sites {input.known_variants} --tmp-dir {resources.tmpdir} -O {output}'
rule ApplyBQSR:
	input:
		bam = lambda x: 'pipeline/mark_dups/{indiv}.downsampled.sorted.mark_dups.{iter}.bam' if 'input_variants' in config.keys() else [],
		reference = config['reference_genome'],
		recal = 'pipeline/BQSR/{indiv}.downsampled.sorted.mark_dups.BQSR.{iter}.txt'
	output:
		temp('pipeline/BQSR/{indiv}.downsampled.sorted.mark_dups.BQSR.{iter}.bam')
	params:
		gatk = config['apps']['gatk']
	resources:
		tmpdir='pipeline/BQSR/'
	shell:
		'{params.gatk} ApplyBQSR -I {input.bam} -R {input.reference} --bqsr-recal-file {input.recal} --tmp-dir {resources.tmpdir} -O {output}'



rule CallVariants:
    input:
        reference = config['reference_genome'],
        input_bam = lambda x: 'pipeline/BQSR/{indiv}.downsampled.sorted.mark_dups.BQSR.{iter}.bam' if 'input_variants' in config.keys() else 'pipeline/mark_dups/{chromosome}_{indiv}.downsampled.sorted.mark_dups.{iter}.bam'
    output:
        output_vcf = temp('pipeline/call_variants/{chromosome}_{indiv}.downsampled.sorted.mark_dups.{iter}.g.vcf'),
	reassembled_bam = 'pipeline/call_variants/{chromosome}_{indiv}.downsampled.sorted.mark_dups.{iter}.reassembled.bam'
    params:
        gatk = config["apps"]["gatk"]
    resources:
    	tmpdir='pipeline/call_variants/',
	runtime='1d',
	mem_mb=get_mem_mb
    shell:
        '{params.gatk} --java-options "-Xmx{resources.mem_mb}m" HaplotypeCaller -R {input.reference}  -I {input.input_bam}  -ERC GVCF --minimum-mapping-quality 40 --max-reads-per-alignment-start 0 --pcr-indel-model NONE -O {output.output_vcf} -L {wildcards.chromosome} --tmp-dir {resources.tmpdir} -bamout {output.reassembled_bam}'

rule GenotypeVariants:
    input:
        p1_vcf='pipeline/call_variants/{chromosome}_parent_1.downsampled.sorted.mark_dups.{iter}.g.vcf',
	p2_vcf='pipeline/call_variants/{chromosome}_parent_2.downsampled.sorted.mark_dups.{iter}.g.vcf',
	ch_vcf='pipeline/call_variants/{chromosome}_child.downsampled.sorted.mark_dups.{iter}.g.vcf'
    output:
        'pipeline/genotype_variants/{chromosome}_trio.downsampled.sorted.mark_dups.{iter}.vcf'
    params:
        reference = config['reference_genome'],
        gatk = config["apps"]["gatk"]
    resources:
    	tmpdir='pipeline/genotype_variants/',
	runtime='1d'
    shell:
        'apptainer exec --bind \"/scratch/mmilhave/\" /scratch/mmilhave/DNM-Coverage_RERUN/images/gatk3_3.7.sif gatk \
-T GenotypeGVCFs \
-R /scratch/mmilhave/refs/panTro6.fa \
--variant {input.p1_vcf} \
--variant {input.p2_vcf} \
--variant {input.ch_vcf} \
-A StrandOddsRatio \
-A MappingQualityRankSumTest \
-A ReadPosRankSumTest \
-A QualByDepth \
-A FisherStrand \
--out {output}'

rule IsolateMVs:
    input:
        'pipeline/genotype_variants/{chromosome}_trio.downsampled.sorted.mark_dups.{iter}.vcf'
    output:
        'pipeline/MV/{chromosome}_trio.downsampled.sorted.mark_dups.MV.{iter}.vcf'
    params:
        gatk = config['apps']['gatk']
    resources:
        runtime='15m' 
    shell:
        '{params.gatk} SelectVariants -V {input} --restrict-alleles-to BIALLELIC --select-type-to-include SNP --exclude-filtered true --exclude-non-variants true -select \'AN==6\' --select \'vc.getGenotype("parent_1").isHomRef()\' --select \'vc.getGenotype("parent_2").isHomRef()\' --select \'vc.getGenotype("child").isHet()\' -O {output}'

rule MergeMVVCFs:
    input:
        expand('pipeline/MV/{chromosome}_trio.downsampled.sorted.mark_dups.MV.{{iter}}.vcf', chromosome=config['chroms'].keys())
    output:
        'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.{iter}.vcf'
    params:
        formatted = format_input,
        gatk = config["apps"]["gatk"]
    resources:
        runtime='15m'
    run:
        shell('{params.gatk} MergeVcfs {params.formatted} -O {output}')

rule get_MV_muts:
	input:
		vcf = 'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.{iter}.vcf',
		muts = 'pipeline/mutations/input_mutations.{iter}.phased.reformatted.vcf'
	output:
		'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.{iter}.mutations.vcf'
	params:
		gatk = config["apps"]["gatk"]
	resources:
		runtime='15m'
	shell:
		'{params.gatk} SelectVariants -V {input.vcf} --concordance {input.muts} -O {output}'

rule get_muts_table:
	input:
		'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.{iter}.mutations.vcf'
	output:
		'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.{iter}.mutations.table'
	params:
		gatk=config["apps"]["gatk"]
	resources:
		runtime='15m'
	shell:
		'{params.gatk} VariantsToTable -V {input} -F CHROM -F POS -F QUAL -F BaseQRankSum -F MQRankSum -F ReadPosRankSum -F SOR -F FS -F QD -GF AD -GF DP -GF GQ -O {output}'

rule get_MV_vars:
	input:
		vcf = 'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.{iter}.vcf',
		poly_vcf = config['input_variants']
	output:
		'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.{iter}.polymorphisms.vcf'
	params:
		gatk = config["apps"]["gatk"]
	resources:
		runtime='15m'
	shell:
		'{params.gatk} SelectVariants -V {input.vcf} --concordance {input.poly_vcf} -O {output}'

rule get_vars_table:
	input:
		'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.{iter}.polymorphisms.vcf'
	output:
		'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.{iter}.polymorphisms.table'
	params:
		gatk = config["apps"]["gatk"]
	resources:
		runtime='15m'
	shell:
		'{params.gatk} VariantsToTable -V {input} -F CHROM -F POS -F QUAL -F BaseQRankSum -F MQRankSum -F ReadPosRankSum -F SOR -F FS -F QD -GF AD -GF DP -GF GQ -O {output}'

rule get_MV_errors:
	input:
		vcf = 'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.{iter}.vcf',
                xvars = lambda x: 'pipeline/mutations/polymorphisms.{iter}.phased.vcf' if 'input_variants' in config.keys() else 0,
		muts = 'pipeline/mutations/input_mutations.{iter}.phased.reformatted.vcf'
	output:
		'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.{iter}.errors.vcf'
	params:
		gatk = config["apps"]["gatk"]
	resources:
		runtime='15m'
	shell:
		'{params.gatk} SelectVariants -V {input.vcf} --discordance {input.muts} -O pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.{wildcards.iter}.nomuts.vcf && {params.gatk} SelectVariants -V pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.{wildcards.iter}.nomuts.vcf --discordance {input.xvars} -O {output}'

rule get_errors_table:
	input:
		'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.{iter}.errors.vcf'
	output:
		'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.{iter}.errors.table'
	params:
		gatk = config["apps"]["gatk"]
	resources:
		runtime='15m'
	shell:
		'{params.gatk} VariantsToTable -V {input} -F CHROM -F POS -F QUAL -F BaseQRankSum -F MQRankSum -F ReadPosRankSum -F SOR -F FS -F QD -GF AD -GF DP -GF GQ -O {output}'

rule get_recommended_stats:
	input:
		mut_table = expand('pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.{iter}.mutations.table',iter=np.arange(1,int(config['num_iterations'])+1)),
		vars_table = expand('pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.{iter}.polymorphisms.table',iter=np.arange(1,int(config['num_iterations'])+1)),
		err_table = expand('pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.{iter}.errors.table',iter=np.arange(1,int(config['num_iterations'])+1))
	output:
		'pipeline/run_outputs/run.{iter}.statistics.txt'
		
	params:
		work_dir = config['working_directory'],
		snakedir=SNAKEDIR,
		num_iterations=config["num_iterations"]
	shell:
		'python {params.snakedir}/scripts/get_recommended_stats.py -n {wildcards.iter}  -o {output}'

rule super_aggregate:
	input:
		expand('pipeline/run_outputs/run{iter}.statistics.txt',iter = range(1,int(config['num_iterations'])+1))
		#get_super_aggregate_input()
	output:
		config['output_file']
	resources:
		runtime='15m'
	params:
		snakedir=SNAKEDIR,
		num_runs=config['num_iterations']
	shell:
		'python {params.snakedir}/scripts/super_aggregate.py -n {params.num_runs} -o {output}'
