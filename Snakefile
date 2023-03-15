#~/.conda/envs/snakemake/bin/python
import os
cwd = os.path.abspath(os.path.dirname(__file__))
SNAKEDIR = config["snakemake_dir"]
workdir: config["working_directory"]
num_chroms = len(config["chroms"].keys())
read_cutoff = 5e6
#TODO Clean up shell commands


import pprint
import numpy
from numpy import random
pp = pprint.PrettyPrinter(indent=4)

wildcard_constraints:
    filter_threshold="-?\d+\.\d+",
    indiv="|".join(list(config['names'].keys())).replace("_", "\_"), # this commands gets all possible value of names from config file and turns it into a regex
    chromosome="|".join(list(config['chroms'].keys())).replace("_", "\_"),
    num="0|1"

def get_all_input():
    target_list = ['pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.mutations.vcf']
    filters = get_filters()
    for filter in filters:
        target_list += expand('pipeline/filters/{filter}/all_chr_trio.downsampled.sorted.mark_dups.MV.{filter}.{filter_threshold}.mutations.vcf', filter_threshold = numpy.round(numpy.arange(float(config["filters"][filter]["min"]),float(config["filters"][filter]["max"])+float(config["filters"][filter]["step"])/2,float(config["filters"][filter]["step"])),2),filter=filter),
        target_list += expand('pipeline/filters/{filter}/all_chr_trio.downsampled.sorted.mark_dups.MV.{filter}.{filter_threshold}.vcf', filter_threshold = numpy.round(numpy.arange(float(config["filters"][filter]["min"]),float(config["filters"][filter]["max"])+float(config["filters"][filter]["step"])/2,float(config["filters"][filter]["step"])),2),filter=filter),
        if 'input_variants' in config.keys():
            target_list += expand('pipeline/filters/{filter}/all_chr_trio.downsampled.sorted.mark_dups.MV.{filter}.{filter_threshold}.polymorphisms.vcf', filter_threshold = numpy.round(numpy.arange(float(config["filters"][filter]["min"]),float(config["filters"][filter]["max"])+float(config["filters"][filter]["step"])/2,float(config["filters"][filter]["step"])),2),filter=filter)
    if 'input_variants' in config.keys():
        target_list.append('pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.polymorphisms.vcf')
    return target_list

def get_filters():
    filters = config['filters']
    return filters.keys()

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

def get_par1_threshold(wildcards):
    return numpy.round(float(wildcards.filter_threshold) * float(config["coverages"]["parent_1"]),2)

def get_par2_threshold(wildcards):
    return numpy.round(float(wildcards.filter_threshold) * float(config["coverages"]["parent_2"]),2)

def get_chd_threshold(wildcards):
    return numpy.round(float(wildcards.filter_threshold) * float(config["coverages"]["child"]),2)

def format_input(wildcards, input):
    return [f'-I {a}' for a in input]

def format_genotype(wildcards,input):
    return [f'--variant {a}' for a in input]

def get_read_count(wildcards): 
    read_length = float(config["read_length"])
    chrom_length = float(config["chroms"][wildcards.chromosome])
    coverage = float(config["coverages"][wildcards.indiv])
    return int(round(((chrom_length*coverage)/read_length)*1.2)/2) 

def get_run_count(wildcards):
    read_count = get_read_count(wildcards)
    to_make = expand('pipeline/sorted_bams/{chrom}_{indiv}_{run}.downsampled.sorted.bam', chrom=wildcards.chromosome, indiv=wildcards.indiv, run=numpy.arange(1, int((read_count / read_cutoff)+1 +1)))
    return to_make

#def get_seeds(wildcards):
#	return numpy.random.randint(1,10000,1)[0]

def get_reads_per_run(wildcards):
    read_count = get_read_count(wildcards)
    num_runs = int((read_count / read_cutoff)+1)
    runs = int(read_count/num_runs)
    return runs

rule all:
    input:
        config['output_file']

rule SplitFasta:
    input:
        fa = config['reference_genome']
    output:
        fa = 'pipeline/split_ref/{chromosome}.fa'
    params:
        samtools = config["apps"]["samtools"]
    shell:
        '{params.samtools} faidx {input.fa} {wildcards.chromosome} > {output.fa}'

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
		'{params.python} {params.snakedir}/scripts/rename_samples.py {input.input_file}  \"{params.p1},{params.p2},{params.ch}\" > {output}'

rule Mutate:                                                                    
    input:                                                                      
        fa = config['reference_genome'],
        vcf = get_Mutate_input if 'input_variants' in config.keys() else [] 
    output:                                                                     
        'pipeline/mutations/all_mutations.vcf'                                  
    params:                                                                     
        mu  = config['mu'],                                                     
        python = config["apps"]["python"],                                      
        chroms = ",".join(list(config["chroms"].keys())),                       
        vcfmerge = format_vcfmerge, 
        snakedir = SNAKEDIR,
	gatk = config["apps"]["gatk"]                                                     
    resources:
        mem_mb = 8000
    run:                                                                      
        shell('{params.python} {params.snakedir}/scripts/mutator.py -i {input.fa} {params.vcfmerge} -u {params.mu} -p \"parent_1,parent_2,child\" -c \"{params.chroms}\" -o {output}')
	shell('{params.gatk} IndexFeatureFile -I {output}')

rule FilterMutations:
    input:
        'pipeline/mutations/all_mutations.vcf'
    output:
        'pipeline/mutations/mutations.vcf'
    params:
        gatk = config["apps"]["gatk"]
    shell:
        '{params.gatk} SelectVariants -V {input} --select "MT==1" -O {output}'

rule FilterPolymorphisms:
    input:
        'pipeline/mutations/all_mutations.vcf'
    params:
        gatk = config["apps"]["gatk"]
    output:
        'pipeline/mutations/polymorphisms.vcf'
    shell:
        '{params.gatk} SelectVariants -V {input} --select "MT==0" -O {output}'

rule SplitInputMutations:
	input:
		input_file = 'pipeline/mutations/all_mutations.vcf',
	output:
		'pipeline/mutations/{chromosome}.mutations.vcf'
	params:
		gatk = config["apps"]["gatk"]
	shell:
		'{params.gatk} SelectVariants -V {input.input_file} -L {wildcards.chromosome} -O {output}'

rule RemoveMutationContigs:
    input:
        'pipeline/mutations/{chromosome}.mutations.vcf'
    output:
        'pipeline/mutations/{chromosome}.mutations.one_contig.vcf'
    params:
        python = config["apps"]["python"],
        snakedir = SNAKEDIR
    shell:
        '{params.python} {params.snakedir}/scripts/remove_contigs.py -i {input} -c {wildcards.chromosome} -o {output}'

rule splitVCFInd:
    input:
        'pipeline/mutations/{chromosome}.mutations.one_contig.vcf'
    output:
        'pipeline/individual_variants/{chromosome}_{indiv}.merged_variants.vcf'
    params:
        gatk = config["apps"]["gatk"]
    shell:
        '{params.gatk} IndexFeatureFile -I {input} && {params.gatk} SelectVariants -V {input} --sample-name {wildcards.indiv} -O {output}'

rule SimulateReads:
    input:
        fa = 'pipeline/split_ref/{chromosome}.fa',
        vcf = 'pipeline/individual_variants/{chromosome}_{indiv}.merged_variants.vcf' if 'input_variants' in config.keys() else [] 
    output:
        r1 = 'pipeline/reads/{chromosome}_{indiv}_{run}_R1.fq', 
        r2 = 'pipeline/reads/{chromosome}_{indiv}_{run}_R2.fq',
    params:
        read_length = config['read_length'],
        read_fragmean = config['read_fragment'],
        mason = config["apps"]["mason_simulator"],
        read_count = get_reads_per_run,
        input_vars = format_simulate_variants, 
    threads:
        int(config['num_cores'])
    shell:
        '{params.mason} --read-name-prefix "simulated_{wildcards.run}" --seed $RANDOM --num-threads {threads} -ir {input.fa} {params.input_vars} --fragment-mean-size {params.read_fragmean} --illumina-read-length {params.read_length} -n {params.read_count} -o {output.r1} -or {output.r2}'

#rule RenameReads:
#	input:
#		r1 = 'pipeline/reads/{chromosome}_{indiv}_{run}_R1.fq',
#		r2 = 'pipeline/reads/{chromosome}_{indiv}_{run}_R2.fq'
#	output:
#		r1 = 'pipeline/reads/{chromosome}_{indiv}_{run}_renamed.R1.fq',
#		r2 = 'pipeline/reads/{chromosome}_{indiv}_{run}_renamed.R2.fq'
#	shell:
#		"sed 's/simulated/simulated_{wildcards.run}/g' {input.r1} > {output.r1} && sed 's/simulated/simulated_{wildcards.run}/g' {input.r2} > {output.r2}"

rule MapReadsDownsampleBam: 
	input:
		reference = config['reference_genome'],
		fq1 = 'pipeline/reads/{chromosome}_{indiv}_{run}_R1.fq', 
		fq2 = 'pipeline/reads/{chromosome}_{indiv}_{run}_R2.fq'
	threads:
		int(config['num_cores'])
	output:
		'pipeline/sorted_bams/{chromosome}_{indiv}_{run}.downsampled.bam'
	params:
		bwa = config["apps"]["bwa"],
		samtools = config['apps']['samtools']
	resources:
		mem_mb = 2000 * int(config['num_cores'])
	shell:
		'{params.bwa} mem -t {threads} -R \"@RG\\tID:{wildcards.chromosome}_{wildcards.run}\\tSM:{wildcards.indiv}\\tPL:ILLUMINA\" {input.reference}  {input.fq1} {input.fq2} | {params.samtools} view -@ {threads} -h -b -s 0.83333333 - > {output}'

rule SortBam:
	input:
		'pipeline/sorted_bams/{chromosome}_{indiv}_{run}.downsampled.bam'
	output:
		'pipeline/sorted_bams/{chromosome}_{indiv}_{run}.downsampled.sorted.bam'
	params:
		samtools = config['apps']['samtools']
	resources:
		mem_mb = 2000 * int(config['num_cores'])
	shell:
		'{params.samtools} sort -@ {threads} {input} > {output}'

rule MergeBAMs:
    input:
        get_run_count
    output:
        'pipeline/sorted_bams/{chromosome}_{indiv}.downsampled.sorted.bam'
    params:
        formatted = format_input,
        gatk = config["apps"]["gatk"]
    resources:
        mem_mb = 8000
    shell:
        '{params.gatk} MergeSamFiles {params.formatted} -O {output}'


rule MarkDuplicates: 
    input:
        'pipeline/sorted_bams/{chromosome}_{indiv}.downsampled.sorted.bam'
    output:
        output_bam = 'pipeline/mark_dups/{chromosome}_{indiv}.downsampled.sorted.mark_dups.bam',
        output_metrics = 'pipeline/mark_dups/{chromosome}_{indiv}.downsampled.sorted.mark_dups.metrics.txt'
    params:
        gatk = config['apps']['gatk']
    shell:
        '{params.gatk} MarkDuplicates -I {input} -O {output.output_bam} -M {output.output_metrics}  --CREATE_INDEX true'

rule BQSR:
	input:
		bam = lambda x: 'pipeline/mark_dups/{chromosome}_{indiv}.downsampled.sorted.mark_dups.bam' if 'input_variants' in config.keys() else [],
		known_variants = lambda x: config['input_variants'] if 'input_variants' in config.keys() else [],
		reference = config['reference_genome']
	output:
		'pipeline/BQSR/{chromosome}_{indiv}.downsampled.sorted.mark_dups.BQSR.txt'
	params:
		gatk = config['apps']['gatk']
	shell:
		'{params.gatk} BaseRecalibrator -R {input.reference} -I {input.bam} --known-sites {input.known_variants} -O {output}'
rule ApplyBQSR:
	input:
		bam = lambda x: 'pipeline/mark_dups/{chromosome}_{indiv}.downsampled.sorted.mark_dups.bam' if 'input_variants' in config.keys() else [],
		reference = config['reference_genome'],
		recal = 'pipeline/BQSR/{chromosome}_{indiv}.downsampled.sorted.mark_dups.BQSR.txt'
	output:
		'pipeline/BQSR/{chromosome}_{indiv}.downsampled.sorted.mark_dups.BQSR.bam'
	params:
		gatk = config['apps']['gatk']
	shell:
		'{params.gatk} ApplyBQSR -I {input.bam} -R {input.reference} --bqsr-recal-file {input.recal} -O {output}'



rule CallVariants:
    input:
        reference = config['reference_genome'],
        input_bam = lambda x: 'pipeline/BQSR/{chromosome}_{indiv}.downsampled.sorted.mark_dups.BQSR.bam' if 'input_variants' in config.keys() else 'pipeline/mark_dups/{chromosome}_{indiv}.downsampled.sorted.mark_dups.bam'
    output:
        output_vcf = 'pipeline/call_variants/{chromosome}_{indiv}.downsampled.sorted.mark_dups.g.vcf',
    params:
        gatk = config["apps"]["gatk"]
    shell:
        '{params.gatk} HaplotypeCaller -R {input.reference}  -I {input.input_bam}  -ERC BP_RESOLUTION --minimum-mapping-quality 40 --max-reads-per-alignment-start 0 --pcr-indel-model NONE -O {output.output_vcf} -L {wildcards.chromosome}'

rule MergeCalls:
    input:
        p1 = 'pipeline/call_variants/{chromosome}_parent_1.downsampled.sorted.mark_dups.g.vcf',
        p2 = 'pipeline/call_variants/{chromosome}_parent_2.downsampled.sorted.mark_dups.g.vcf',
        ch = 'pipeline/call_variants/{chromosome}_child.downsampled.sorted.mark_dups.g.vcf' 
    output:
        outdir = directory('pipeline/GBDImport/{chromosome}')
    params:
        gatk = config["apps"]["gatk"]
    shell:
        '{params.gatk} GenomicsDBImport -V {input.p1} -V {input.p2} -V {input.ch} --genomicsdb-workspace-path {output.outdir} -L {wildcards.chromosome}'

rule GenotypeVariants:
    input:
        indir='pipeline/GBDImport/{chromosome}'
    output:
        'pipeline/genotype_variants/{chromosome}_trio.downsampled.sorted.mark_dups.vcf'
        
    params:
        reference = config['reference_genome'],
        gatk = config["apps"]["gatk"]
    shell:
        '{params.gatk} GenotypeGVCFs -R {params.reference} -V gendb://{input.indir} -A StrandOddsRatio -A MappingQualityRankSumTest -A ReadPosRankSumTest -A QualByDepth -A FisherStrand -O {output}'

rule IsolateMVs:
    input:
        'pipeline/genotype_variants/{chromosome}_trio.downsampled.sorted.mark_dups.vcf'
    output:
        'pipeline/MV/{chromosome}_trio.downsampled.sorted.mark_dups.MV.vcf'
    params:
        gatk = config['apps']['gatk']
        
    shell:
        '{params.gatk} SelectVariants -V {input} --restrict-alleles-to BIALLELIC --select-type-to-include SNP --exclude-filtered true --exclude-non-variants true -select \'AN==6\' --select \'vc.getGenotype("parent_1").isHomRef()\' --select \'vc.getGenotype("parent_2").isHomRef()\' --select \'vc.getGenotype(\"child\").isHet()\' -O {output}'

rule MergeMVVCFs:
    input:
        expand('pipeline/MV/{chromosome}_trio.downsampled.sorted.mark_dups.MV.vcf', chromosome=config['chroms'].keys())
    output:
        'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.vcf'
    params:
        formatted = format_input,
        gatk = config["apps"]["gatk"]
    run:
        shell('{params.gatk} MergeVcfs {params.formatted} -O {output}')

rule DPGT_filter:
    input:
        'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.vcf'
    output:
        'pipeline/filters/DPGT/all_chr_trio.downsampled.sorted.mark_dups.MV.DPGT.{filter_threshold}.vcf'
    params:
        parent1_threshold = get_par1_threshold,
        parent2_threshold = get_par2_threshold,
        child_threshold = get_chd_threshold,
        gatk = config["apps"]["gatk"]
    shell:
        '{params.gatk} SelectVariants -V {input} --select \'vc.getGenotype(\"parent_1\").getDP() >= {params.parent1_threshold}\' --select \'vc.getGenotype(\"parent_2\").getDP() >= {params.parent2_threshold}\' --select \'vc.getGenotype(\"child\").getDP() >= {params.child_threshold}\' --exclude-filtered true -O {output}'

rule DPLT_filter:
    input:
        'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.vcf'
    output:
        'pipeline/filters/DPLT/all_chr_trio.downsampled.sorted.mark_dups.MV.DPLT.{filter_threshold}.vcf'
    params:        
        parent1_threshold = get_par1_threshold,
        parent2_threshold = get_par2_threshold,
        child_threshold = get_chd_threshold,
        gatk = config["apps"]["gatk"]
    shell:                                                                      
        '{params.gatk} SelectVariants -V {input} --select \'vc.getGenotype(\"parent_1\").getDP() <= {params.parent1_threshold}\' --select \'vc.getGenotype(\"parent_2\").getDP() <= {params.parent2_threshold}\' --select \'vc.getGenotype(\"child\").getDP() <= {params.child_threshold}\' --exclude-filtered true -O {output} '
 
rule ABGT_filter:
    input:
        'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.vcf'
    output:
        'pipeline/filters/ABGT/all_chr_trio.downsampled.sorted.mark_dups.MV.ABGT.{filter_threshold}.vcf'
    params:
        gatk = config["apps"]["gatk"]
    shell:
        '{params.gatk} SelectVariants -V {input} --select \'(vc.getGenotype(\"child\").getAD().1*1.0 / vc.getGenotype(\"child\").getDP()*1.0 ) >= {wildcards.filter_threshold} \' --exclude-filtered true -O {output}'

rule ABLT_filter:
    input:
        'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.vcf'
    output:
        'pipeline/filters/ABLT/all_chr_trio.downsampled.sorted.mark_dups.MV.ABLT.{filter_threshold}.vcf'
    params:
        gatk = config["apps"]["gatk"]
    shell:
        '{params.gatk} SelectVariants -V {input} --select \'(vc.getGenotype(\"child\").getAD().1*1.0 / vc.getGenotype(\"child\").getDP()*1.0 ) <= {wildcards.filter_threshold} \' --exclude-filtered true -O {output}'

rule QUAL_filter:
    input:
        'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.vcf'
    output:
        'pipeline/filters/QUAL/all_chr_trio.downsampled.sorted.mark_dups.MV.QUAL.{filter_threshold}.vcf'
    params:
        gatk = config["apps"]["gatk"]
    shell:
        '{params.gatk} SelectVariants -V {input} --select \"QUAL >= {wildcards.filter_threshold}\"  --exclude-filtered true -O {output}'

rule GQ_filter:
    input:
        'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.vcf'
    output:
        'pipeline/filters/GQ/all_chr_trio.downsampled.sorted.mark_dups.MV.GQ.{filter_threshold}.vcf'
    params:
        gatk = config["apps"]["gatk"]
    shell:
        '{params.gatk} SelectVariants -V {input} --select \'vc.getGenotype(\"parent_1\").getGQ() >= {wildcards.filter_threshold}\' --select \'vc.getGenotype(\"parent_2\").getGQ() >= {wildcards.filter_threshold}\' --select \'vc.getGenotype(\"child\").getGQ() >= {wildcards.filter_threshold}\' --exclude-filtered true -O {output} '

rule AD_filter:
    input:
        'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.vcf'
    output:
        'pipeline/filters/AD/all_chr_trio.downsampled.sorted.mark_dups.MV.AD.{filter_threshold}.vcf'
    params:
        gatk = config["apps"]["gatk"]
    shell:
        '{params.gatk} SelectVariants -V {input} --select \'vc.getGenotype(\"parent_1\").getAD().1 <= {wildcards.filter_threshold}\' --select \'vc.getGenotype(\"parent_2\").getAD().1 <= {wildcards.filter_threshold} \' --exclude-filtered true -O {output} '

#rule GATKBP_filter:
#    input:
#        'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.vcf'
#    output:
#        gatkbp_output = 'pipeline/filters/GATKBP/all_chr_trio.downsampled.sorted.mark_dups.MV.GATKBP.vcf',
#        gatkbp_pass_output = 'pipeline/filters/GATKBP/all_chr_trio.downsampled.sorted.mark_dups.MV.GATKBP.1.0.vcf'
#    params:
#        gatk = config["apps"]["gatk"]
#    run:
#        shell('{params.gatk} VariantFiltration -V {input} -filter \"QD < 2.0\" --filter-name \"QD2\" -filter \"SOR > 3.0\" --filter-name \"SOR3\" -filter \"FS > 60.0\" --filter-name \"FS60\" -filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" -filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\" -O {output.gatkbp_output}')
#        shell('{params.gatk} SelectVariants -V {output.gatkbp_output} --exclude-filtered true -O {output.gatkbp_pass_output}')

rule QD_filter:
	input:
		'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.vcf'
	output:
		'pipeline/filters/QD/all_chr_trio.downsampled.sorted.mark_dups.MV.QD.{filter_threshold}.vcf'
	params:
		gatk = config["apps"]["gatk"]
	shell:
		'{params.gatk} SelectVariants -V {input} --select "QD >= {wildcards.filter_threshold}" -O {output}'

rule FS_filter:
	input:
		'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.vcf'
	output:
		'pipeline/filters/FS/all_chr_trio.downsampled.sorted.mark_dups.MV.FS.{filter_threshold}.vcf'
	params:
		gatk = config["apps"]["gatk"]
	shell:
		'{params.gatk} SelectVariants -V {input} --select "FS <= {wildcards.filter_threshold}" -O {output}'

rule SOR_filter:
	input:
		'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.vcf'
	output:
		'pipeline/filters/SOR/all_chr_trio.downsampled.sorted.mark_dups.MV.SOR.{filter_threshold}.vcf'
	params:
		gatk = config["apps"]["gatk"]
	shell:
		'{params.gatk} SelectVariants -V {input} --select "SOR <= {wildcards.filter_threshold}" -O {output}'

rule MQRankSum_GT_filter:
	input:
		'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.vcf'
	output:
		'pipeline/filters/MQRankSumGT/all_chr_trio.downsampled.sorted.mark_dups.MV.MQRankSumGT.{filter_threshold}.vcf'	
	params:
		gatk = config["apps"]["gatk"]
	shell:
		'{params.gatk} SelectVariants -V {input} --select "MQRankSum >= {wildcards.filter_threshold}" -O {output}'

rule MQRankSum_LT_filter:
	input:
		'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.vcf'
	output:
		'pipeline/filters/MQRankSumLT/all_chr_trio.downsampled.sorted.mark_dups.MV.MQRankSumLT.{filter_threshold}.vcf'
	params:
		gatk = config["apps"]["gatk"]
	shell:
		'{params.gatk} SelectVariants -V {input} --select "MQRankSum <= {wildcards.filter_threshold}" -O {output}'

rule ReadPosRankSum_GT_filter:
	input:
		'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.vcf'
	output:
		'pipeline/filters/ReadPosRankSumGT/all_chr_trio.downsampled.sorted.mark_dups.MV.ReadPosRankSumGT.{filter_threshold}.vcf'
	params:
		gatk = config["apps"]["gatk"]
	shell:
		'{params.gatk} SelectVariants -V {input} --select "ReadPosRankSum >= {wildcards.filter_threshold}" -O {output}'

rule ReadPosRankSum_LT_filter:
	input:
		'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.vcf'
	output:
		'pipeline/filters/ReadPosRankSumLT/all_chr_trio.downsampled.sorted.mark_dups.MV.ReadPosRankSumLT.{filter_threshold}.vcf'
	params:
		gatk = config["apps"]["gatk"]
	shell:
		'{params.gatk} SelectVariants -V {input} --select "ReadPosRankSum <= {wildcards.filter_threshold}" -O {output}'

rule get_MV_muts:
	input:
		vcf = 'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.vcf',
		muts = 'pipeline/mutations/mutations.vcf'
	output:
		'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.mutations.vcf'
	params:
		gatk = config["apps"]["gatk"]
	shell:
		'{params.gatk} SelectVariants -V {input.vcf} --concordance {input.muts} -O {output}'

rule get_MV_vars:
	input:
		vcf = 'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.vcf',
		xvars = lambda x: 'pipeline/mutations/polymorphisms.vcf' if 'input_variants' in config.keys() else 0
	output:
		'pipeline/MV/all_chr_trio.downsampled.sorted.mark_dups.MV.polymorphisms.vcf'
	params:
		gatk = config["apps"]["gatk"]
	shell:
		'{params.gatk} SelectVariants -V {input.vcf} --concordance {input.xvars} -O {output}'
rule get_muts:
    input:
        vcf = 'pipeline/filters/{filter}/all_chr_trio.downsampled.sorted.mark_dups.MV.{filter}.{filter_threshold}.vcf',
        muts = 'pipeline/mutations/mutations.vcf'
    output:
        'pipeline/filters/{filter}/all_chr_trio.downsampled.sorted.mark_dups.MV.{filter}.{filter_threshold}.mutations.vcf'
    params:
        gatk = config["apps"]["gatk"]
    shell:
        '{params.gatk} SelectVariants -V {input.vcf} --concordance {input.muts} -O {output}'

rule get_vars:
    input:
        vcf = 'pipeline/filters/{filter}/all_chr_trio.downsampled.sorted.mark_dups.MV.{filter}.{filter_threshold}.vcf',
        xvars = lambda x: 'pipeline/mutations/polymorphisms.vcf' if 'input_variants' in config.keys() else 0
    output:
        'pipeline/filters/{filter}/all_chr_trio.downsampled.sorted.mark_dups.MV.{filter}.{filter_threshold}.polymorphisms.vcf'
    params:
        gatk = config["apps"]["gatk"]
    shell:
        '{params.gatk} SelectVariants -V {input.vcf} --concordance {input.xvars} -O {output}'

rule aggregate:
    input:
        get_all_input()
    output:
        config['output_file']
    params:
        vcfmerge = lambda x: 1 if 'input_variants' in config.keys() else 0,
        snakedir=SNAKEDIR,
        config = config["working_directory"]+"/config/config.json", #["filter_file"],
        work_dir= config["working_directory"]
    shell:
        'python {params.snakedir}/scripts/aggregate.py -i {params.config} -o {output} -v {params.vcfmerge} -m pipeline/mutations/multihit_count.txt -w {params.work_dir}'

