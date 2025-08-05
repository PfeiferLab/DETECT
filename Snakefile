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
    indiv="|".join(list(config['names'].keys())).replace("_", "\\_"), 
    chromosome="|".join(list(config['chroms'].keys())).replace("_", "\\_"),
    genome="0|1",
    iter='|'.join([str(x) for x in numpy.arange(1,9999)]),
    sim="|".join([str(x) for x in numpy.arange(1,9999)])

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

def format_sortbam_input(wildcards):
	if wildcards.indiv == 'child':
		return 'pipeline/sorted_bams/child.downsampled.{iter}.bam'
	else:
		return 'pipeline/sorted_bams/{indiv}.{iter}.bam'
def format_genotype(wildcards,input):
    return [f'--variant {a}' for a in input]

def get_read_count(wildcards): 
	read_length = float(config["read_length"])
	chrom_length = float(sum([int(x) for x in config["chroms"].values()]))
	coverage = float(config["coverages"][wildcards.indiv])
	if wildcards.indiv == 'child':
		return int(round(((chrom_length*coverage)/read_length)*1.2)/2) 
	else:
		return int(round(((chrom_length*coverage)/read_length))/2)

def get_reads_per_sim(wildcards):
    read_count = get_read_count(wildcards)
    num_sims = 16
    num_reads = int(read_count/num_sims)
    return num_reads

def get_simreads_runtime(wildcards):
	num_reads = get_reads_per_run(wildcards)
	runtime = str(26.9*(num_reads/1e6) + 237 + 60)+'m'
	return runtime

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
		snakedir = SNAKEDIR
	shell:
		'python {params.snakedir}/scripts/rename_samples.py {input.input_file}  \"{params.p1},{params.p2},{params.ch}\" > {output} && gatk IndexFeatureFile -I {output}'

rule Mutate:
	input:
		fa = config['reference_genome']
	output:
		muts='pipeline/mutations/input_mutations.{iter}.phased.vcf'
	params:
		mu  = config['dnm_file'] if config['dnm_file'] != "False" else config['mu'],
		snakedir = SNAKEDIR
	resources:
		mem_mb = 8000,
		runtime='15m'
	shell:
		'python {params.snakedir}/scripts/mutator.py -i {input.fa} -u {params.mu} -p \"parent_1,parent_2,child\" -n {wildcards.iter} -o {output.muts}; '
	'gatk IndexFeatureFile -I {output.muts}'

rule ReformatMutations:
	input:
		muts = config['dnm_file'] if config['dnm_file'] != 'False' else 'pipeline/mutations/input_mutations.{iter}.phased.vcf'
	output:
		'pipeline/mutations/input_mutations.{iter}.phased.reformatted.vcf'
	params:
		snakedir=config["snakemake_dir"],
	resources:
		runtime='15m'
	shell:
		'python {params.snakedir}/scripts/reformat_vcf.py {input.muts} > {output} && gatk IndexFeatureFile -I {output}'

rule splitVCF:
	input:
		muts='pipeline/input_variants/input_trio.renamed.vcf'
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
		'pipeline/ref/parent_1_{chromosome}:0.fa',
		'pipeline/ref/parent_1_{chromosome}:1.fa',
		'pipeline/ref/parent_2_{chromosome}:0.fa',
		'pipeline/ref/parent_2_{chromosome}:1.fa',
		'pipeline/ref/child_{chromosome}:0.fa',
		'pipeline/ref/child_{chromosome}:1.fa' 
	params:
		ref = config['reference_genome'],
	resources:
		runtime='4d'
	shell:
		'mkdir -p pipeline/ref && cd pipeline/ref && vcf2fasta -f {input.ref} ../mutations/polymorphisms.{wildcards.chromosome}.phased.vcf'

rule MergeFastas:
	input:
		expand('pipeline/ref/{{indiv}}_{chromosome}:{{genome}}.fa',chromosome=config['chroms'].keys())
	output:
		'pipeline/ref/{indiv}_renamed.{genome}.fa'
	shell:
		'cat {input} | sed "s/:{wildcards.genome}//g; s/{wildcards.indiv}_//g" > {output}'

rule SimulateReads:
	input:
		fa = 'pipeline/ref/{indiv}_renamed.{genome}.fa',
	output:
		r1 = temp('pipeline/reads/{indiv}_{sim}_{iter}_{genome}.golden.R1.fq'),
		r2 = temp('pipeline/reads/{indiv}_{sim}_{iter}_{genome}.golden.R2.fq'),
		golden_bam = temp('pipeline/reads/{indiv}_{sim}_{iter}_{genome}.golden.bam')
	params:
		read_length = config['read_length'],
		read_fragmean = config['read_fragment'],
		read_count = get_reads_per_sim,
	retries: 20
	resources:
		mem_mb=get_mem_mb,
		runtime='2d'#get_simreads_runtime
	threads:
		int(config['num_cores'])
	shell:
		'mason_simulator --read-name-prefix "simulated_{wildcards.indiv}_{wildcards.sim}_{wildcards.iter}_{wildcards.genome}" --seed $RANDOM --num-threads 12 -ir {input.fa} --fragment-mean-size {params.read_fragmean} --illumina-read-length {params.read_length} -n {params.read_count} -o {output.r1} -or {output.r2} -oa {output.golden_bam}' 

rule MergeGoldenBam:
	input:
		expand('pipeline/reads/{{indiv}}_{sim}_{{iter}}_{genome}.golden.bam',sim=numpy.arange(1,9,1),genome=[0,1])
	params:
		input_list = format_merge_input,
	output:
		temp('pipeline/reads/{indiv}_{iter}.golden.bam')
	resources:
		runtime='8h'
	threads:
		int(config['num_cores'])
	shell:
		'samtools merge -@ {threads} {output} {params.input_list}'

rule MutateGoldenBam:
	input:
		bam='pipeline/reads/child_{iter}.golden.bam',
		vcf='pipeline/mutations/input_mutations.{iter}.phased.reformatted.vcf'
	output:
		temp('pipeline/reads/child_{iter}.golden.mutated.bam')
	resources:
		runtime='4d'
	shell:
		'jvarkit biostar404363 -p {input.vcf} -o {output} {input.bam}'

rule SamToFastq:
	input:
		format_samtofastq_input
	output:
		fq1=temp('pipeline/reads/{indiv}_{iter}.R1.fq'),
		fq2=temp('pipeline/reads/{indiv}_{iter}.R2.fq')
	shell:
		'gatk SamToFastq -I {input} -F {output.fq1} -F2 {output.fq2}'

rule MapReads: 
	input:
		reference = config['reference_genome'],
		fq1 = 'pipeline/reads/{indiv}_{iter}.R1.fq', 
		fq2 = 'pipeline/reads/{indiv}_{iter}.R2.fq',
	threads: int(config['num_cores'])
	output:
		temp('pipeline/sorted_bams/{indiv}.{iter}.bam')
	resources:
		mem_mb = 2000 * int(config['num_cores']),
		runtime='2d'
	shell:
		'bwa mem -t {threads} -R \"@RG\\tID:{wildcards.indiv}_{wildcards.iter}\\tSM:{wildcards.indiv}\\tPL:ILLUMINA\" {input.reference}  {input.fq1} {input.fq2} | samtools view -b - > {output} '

rule DownsampleBam:
	input:
		'pipeline/sorted_bams/child.{iter}.bam'
	output:
		temp('pipeline/sorted_bams/child.downsampled.{iter}.bam')
	resources:
		mem_mb = 2000 * int(config['num_cores']),
		runtime='4h'
	threads: int(config['num_cores'])
	
	shell:
		'samtools view -@ {threads} -b -s 0.83333333 {input} > {output}'

rule SortBam:
	input:
		format_sortbam_input
	output:
		temp('pipeline/sorted_bams/{indiv}.sorted.{iter}.bam')
	resources:
		mem_mb=lambda wc, input: max(input.size_mb*2,1024),
		runtime='12h'
	shell:
		'gatk SortSam --java-options "-Xmx{resources.mem_mb}m" -I {input} --TMP_DIR pipeline/sorted_bams/ -O {output} -SO coordinate '
	

rule MarkDuplicates:
	input:
    		'pipeline/sorted_bams/{indiv}.sorted.{iter}.bam'
	output:
		output_bam = temp('pipeline/mark_dups/{indiv}.sorted.mark_dups.{iter}.bam'),
		output_bai = temp('pipeline/mark_dups/{indiv}.sorted.mark_dups.{iter}.bai'),
		output_metrics = temp('pipeline/mark_dups/{indiv}.sorted.mark_dups.metrics.{iter}.txt')
	resources:
		runtime='2d',
		mem_mb=122880 #lambda wc, input: input.size_mb*2,
	shell:
		'gatk --java-options "-Xmx{resources.mem_mb}m" MarkDuplicates -I {input} -O {output.output_bam} -M {output.output_metrics}  --TMP_DIR pipeline/mark_dups/ --CREATE_INDEX true'

rule BQSR:
	input:
		bam = lambda x: 'pipeline/mark_dups/{indiv}.sorted.mark_dups.{iter}.bam' if 'input_variants' in config.keys() else [],
		known_variants = lambda x: config['input_variants'] if 'known_variants' == "NONE" else config['known_variants'],
		reference = config['reference_genome']
	output:
		'pipeline/BQSR/{indiv}.sorted.mark_dups.BQSR.{iter}.txt'
	resources:
		tmpdir='pipeline/BQSR/',
		runtime='8h'
	shell:
		'gatk BaseRecalibrator -R {input.reference} -I {input.bam} --known-sites {input.known_variants} --tmp-dir {resources.tmpdir} -O {output}'
rule ApplyBQSR:
	input:
		bam = 'pipeline/mark_dups/{indiv}.sorted.mark_dups.{iter}.bam',
		reference = config['reference_genome'],
		recal = 'pipeline/BQSR/{indiv}.sorted.mark_dups.BQSR.{iter}.txt'
	output:
		output_bam = 'pipeline/BQSR/{indiv}.sorted.mark_dups.BQSR.{iter}.bam',
		output_bai = 'pipeline/BQSR/{indiv}.sorted.mark_dups.BQSR.{iter}.bai'
	resources:
		tmpdir='pipeline/BQSR/',
		runtime='12h'
	shell:
		'gatk ApplyBQSR -I {input.bam} -R {input.reference} --bqsr-recal-file {input.recal} --tmp-dir {resources.tmpdir} -O {output.output_bam} --create-output-bam-index true'



rule CallVariants:
	input:
		reference = config['reference_genome'],
		input_bam = 'pipeline/BQSR/{indiv}.sorted.mark_dups.BQSR.{iter}.bam',
		input_bai = 'pipeline/BQSR/{indiv}.sorted.mark_dups.BQSR.{iter}.bai'
	output:
		output_vcf = 'pipeline/call_variants/{indiv}.sorted.mark_dups.BQSR.{chromosome}.{iter}.g.vcf',
		reassembled_bam = 'pipeline/call_variants/{indiv}.sorted.mark_dups.BQSR.{chromosome}.{iter}.reassembled.bam'
	resources:
		tmpdir='pipeline/call_variants/',
		runtime='8h',
		mem_mb=40960
	shell:
		'gatk --java-options "-Xmx{resources.mem_mb}m" HaplotypeCaller -R {input.reference}  -I {input.input_bam}  -ERC GVCF --minimum-mapping-quality 40 --max-reads-per-alignment-start 0 --pcr-indel-model NONE -O {output.output_vcf} -L {wildcards.chromosome} --tmp-dir {resources.tmpdir} -bamout {output.reassembled_bam}'

rule CombineVariants:
	input:
		reference = config['reference_genome'],
		p1_vcf='pipeline/call_variants/parent_1.sorted.mark_dups.BQSR.{chromosome}.{iter}.g.vcf',
		p2_vcf='pipeline/call_variants/parent_2.sorted.mark_dups.BQSR.{chromosome}.{iter}.g.vcf',
		ch_vcf='pipeline/call_variants/child.sorted.mark_dups.BQSR.{chromosome}.{iter}.g.vcf'
	output:
		output_dir = directory('pipeline/genotype_variants/trio.{chromosome}.{iter}')
	resources:
		runtime='8h',
		mem_mb=40960
	shell:
		'gatk --java-options "-Xmx{resources.mem_mb}m" \
			GenomicsDBImport \
			-R {input.reference} \
			-V {input.p1_vcf} \
			-V {input.p2_vcf} \
			-V {input.ch_vcf} \
			--genomicsdb-workspace-path {output.output_dir} \
			-L {wildcards.chromosome}'

rule GenotypeVariants:
	input:
		input_dir = 'pipeline/genotype_variants/trio.{chromosome}.{iter}',
		reference = config['reference_genome']
	output:
		'pipeline/genotype_variants/trio.{chromosome}.{iter}.vcf'
	resources:
		tmpdir='pipeline/genotype_variants/',
		runtime='4h',
		mem_mb=40960
	shell:
		'gatk --java-options \"-Xmx{resources.mem_mb}m\" GenotypeGVCFs \
-R {input.reference} \
-V gendb://{input.input_dir} \
-A StrandOddsRatio \
-A MappingQualityRankSumTest \
-A ReadPosRankSumTest \
-A QualByDepth \
-A FisherStrand \
-O {output}'

rule IsolateMVs:
    input:
        'pipeline/genotype_variants/trio.{chromosome}.{iter}.vcf'
    output:
        'pipeline/MV/trio.MV.{chromosome}.{iter}.vcf'
    resources:
        runtime='15m' 
    shell:
        'gatk SelectVariants -V {input} --restrict-alleles-to BIALLELIC --select-type-to-include SNP --exclude-filtered true --exclude-non-variants true -select \'AN==6\' --select \'vc.getGenotype("parent_1").isHomRef()\' --select \'vc.getGenotype("parent_2").isHomRef()\' --select \'vc.getGenotype("child").isHet()\' -O {output}'

rule MergeMVVCFs:
    input:
        expand('pipeline/MV/trio.MV.{chromosome}.{{iter}}.vcf', chromosome=config['chroms'].keys())
    output:
        'pipeline/MV/trio.MV.all_chr.{iter}.vcf'
    params:
        formatted = format_input,
    resources:
        runtime='15m'
    run:
        shell('gatk MergeVcfs {params.formatted} -O {output}')

rule get_MV_muts:
	input:
		vcf = 'pipeline/MV/trio.MV.all_chr.{iter}.vcf',
		muts = 'pipeline/mutations/input_mutations.{iter}.phased.reformatted.vcf'
	output:
		'pipeline/MV/trio.MV.all_chr.{iter}.mutations.vcf'
	resources:
		runtime='15m'
	shell:
		'gatk SelectVariants -V {input.vcf} --concordance {input.muts} -O {output}'

rule get_assembly_depths:
	input:
		input_vcf='pipeline/MV/trio.MV.all_chr.{iter}.mutations.vcf',
		all_reassembled = expand('pipeline/call_variants/{indiv}.sorted.mark_dups.BQSR.{chromosome}.{{iter}}.reassembled.bam',chromosome=config['chroms'].keys(),indiv=config['names'].keys()),
		ch_bqsr_bam = 'pipeline/BQSR/child.sorted.mark_dups.BQSR.{iter}.bam'
	output:
		'pipeline/MV/trio.MV.all_chr.{iter}.assembly_depths.txt'
	resources:
		runtime='4h'
	shell:r"""
{{
    echo "CHROM POS p1_reassembled p2_reassembled child_reassembled"
        grep -v "#" {input.input_vcf} | awk '{{print $1"\t"$2}}' | while read chrom pos
    do
        p1_reassembled=$(samtools depth -aa -r ${{chrom}}:${{pos}}-${{pos}} pipeline/call_variants/${{chrom}}_parent_1.sorted.mark_dups.{wildcards.iter}.reassembled.bam | awk '{{print $3}}')
        p2_reassembled=$(samtools depth -aa -r ${{chrom}}:${{pos}}-${{pos}} pipeline/call_variants/${{chrom}}_parent_2.sorted.mark_dups.{wildcards.iter}.reassembled.bam | awk '{{print $3}}')
        ch_bqsr=$(samtools depth -aa -r ${{chrom}}:${{pos}}-${{pos}} {input.ch_bqsr_bam} | awk '{{print $3}}')
        ch_reassembled=$(samtools depth -aa -r ${{chrom}}:${{pos}}-${{pos}} pipeline/call_variants/${{chrom}}_child.sorted.mark_dups.{wildcards.iter}.reassembled.bam | awk '{{print $3}}')

        if [ "$ch_bqsr" -gt "$ch_reassembled" ]; then
            ch_diff=$((ch_bqsr - ch_reassembled))
        else
            ch_diff=$((ch_reassembled - ch_bqsr))
        fi

        echo "$chrom $pos $p1_reassembled $p2_reassembled $ch_diff"
    done
}} > {output}
"""

rule get_muts_table:
	input:
		'pipeline/MV/trio.MV.all_chr.{iter}.mutations.vcf'
	output:
		'pipeline/MV/trio.MV.all_chr.{iter}.mutations.table'
	resources:
		runtime='15m'
	shell:
		'gatk VariantsToTable -V {input} -F CHROM -F POS -F QUAL -F BaseQRankSum -F MQRankSum -F ReadPosRankSum -F SOR -F FS -F QD -GF AD -GF DP -GF GQ -GF PL -O {output}'

rule get_MV_vars:
	input:
		vcf = 'pipeline/MV/trio.MV.all_chr.{iter}.vcf',
		poly_vcf = config['input_variants']
	output:
		'pipeline/MV/trio.MV.all_chr.{iter}.polymorphisms.vcf'
	resources:
		runtime='15m'
	shell:
		'gatk SelectVariants -V {input.vcf} --concordance {input.poly_vcf} -O {output}'

rule get_vars_table:
	input:
		'pipeline/MV/trio.MV.all_chr.{iter}.polymorphisms.vcf'
	output:
		'pipeline/MV/trio.MV.all_chr.{iter}.polymorphisms.table'
	resources:
		runtime='15m'
	shell:
		'gatk VariantsToTable -V {input} -F CHROM -F POS -F QUAL -F BaseQRankSum -F MQRankSum -F ReadPosRankSum -F SOR -F FS -F QD -GF AD -GF DP -GF GQ -GF PL -O {output}'

rule get_MV_errors:
	input:
		vcf = 'pipeline/MV/trio.MV.all_chr.{iter}.vcf',
		muts = 'pipeline/mutations/input_mutations.{iter}.phased.reformatted.vcf',
		xvars = config['input_variants']
	output:
		'pipeline/MV/trio.MV.all_chr.{iter}.errors.vcf'
	resources:
		runtime='15m'
	shell:
		'gatk SelectVariants -V {input.vcf} --discordance {input.muts} -O pipeline/MV/trio.MV.all_chr.{wildcards.iter}.nomuts.vcf && gatk SelectVariants -V pipeline/MV/trio.MV.{wildcards.iter}.nomuts.vcf --discordance {input.xvars} -O {output}'

rule get_errors_table:
	input:
		input_vcf = 'pipeline/MV/trio.MV.all_chr.{iter}.errors.vcf',
	output:
		'pipeline/MV/trio.MV.all_chr.{iter}.errors.table'
	resources:
		runtime='15m'
	shell:
		'gatk VariantsToTable -V {input.input_vcf} -F CHROM -F POS -F QUAL -F BaseQRankSum -F MQRankSum -F ReadPosRankSum -F SOR -F FS -F QD -GF AD -GF DP -GF GQ -GF PL -O {output}'

rule make_output:
	input:
		mutation_list=expand('pipeline/MV/trio.MV.all_chr.{iter}.mutations.table',iter=np.arange(1,int(config['num_iterations'])+1)),
		reassembly_file=expand('pipeline/MV/trio.MV.all_chr.{iter}.assembly_depths.txt',iter=np.arange(1,int(config['num_iterations'])+1))
	output:
		config['output_file']
	resources:
		runtime='15m'
	params:
		snakedir=SNAKEDIR
	shell:
		'python {params.snakedir}/scripts/make_output_file.py -i {input.mutation_list} -r {input.reassembly_file} -o {output}'

