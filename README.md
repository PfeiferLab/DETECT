# DETECT (DNM Extraction Through Empirical Cutoff Thresholds)

## Description:
DETECT is a simulation-based workflow that recommends filter thresholds in direct mutation rate estimation. By populating DNMs in a simulated trio at a specified mutation rate, we can determine the filter that isolates as many DNMs as possible in the simulation, while also limiting the number of False Positives(FPs). DETECT also takes variant datasets as input to provide an accurate estimation of FPs that the workflow then must discern from true DNMs. This has only been tested and is only functional on diploid, sexually reproducing organisms, but plans to allow asexual reproduction are being worked on.

## Setting Up:
### Environment Installation:
`conda env create -n detect-env -f DETECT/detect_env.yml`  
`source activate detect-env`  

**Note:** the environment name(-n) can be whatever name you would like it to be

### Required Software:
Parenthetical versions refer to the versions that DETECT was tested with. Older/newer versions may work as well.
python3 - https://www.python.org/
Mason(v2.0.9) - https://github.com/seqan/seqan/blob/master/apps/mason2/
GATK4 - https://github.com/broadinstitute/gatk/releases
bwa(v0.7.17) - https://github.com/lh3/bwa
samtools(v1.9) - http://www.htslib.org/download/


## Quickstart:
### Required Inputs:
**Reference Genome:**  
* Reference genome to be used in your real data workflow, where the simulated read data will come from. Must have a dictionary file (GATK CreateSequenceDictionary) and be bwa indexed (bwa index ) ex. reference.fa  
**Mutation Rate:**  
* Estimated mutation rate  of the dataset in question. Can be in scientific notation. ex. 1e-8  
**Filter File:**  
* List of filters you wish to apply to the simulated dataset (format described below). ex. filter_file.txt  
**Read Length:**  
* Length of the reads used in the real dataset. ex. 100  
**Coverage:**  
* Comma-delimited string of the coverages of the sire,dam,offspring. ex. “30,40,50”  
**Output File:**  
* Output file name of consolidated filter recommendations.  

###Optional Inputs:  
**Input Variants:**  
* VCF file containing variants to be used as False Positives. Must have either --trio or --population specified. If --trio, --pedigree required, and only the trio can be in the VCF. Must be indexed (e.g. GATK IndexFeatureFile). If --population, DETECT will “create” an offspring from two random individuals’ haplotypes.  
**Pedigree:**  
* Comma delimited string of the names of sire, dam, and offspring in the VCF (ex. “dad,mom,junior”)  
**Fragment Length:**  
* Mean length of the fragment size distribution of the real data. Default: 300 (Mason Default)  
**Fragment Length Standard Deviation:**  
* Standard deviation of the fragment size distribution of the real data. Default: 30 (Mason Default)  
**Chromosome list:**  
* File of chromosome names to be simulated, one per line. Default: All contigs  
**Path list:**  
* Space-delimited text file that shows the native paths of each of the required. DETECT will use the native command by default (e.g. “gatk” or “samtools”).  
Example of path list:  
```
samtools /packages/apps/spack/18/opt/spack/gcc-12.1.0/samtools-1.9-arv/bin/samtools
gatk /packages/apps/spack/18/opt/spack/gcc-11.2.0/gatk-4.2.6.1-3ds/bin/gatk  
mason_simulator ~/mason/bin/mason_simulator
```  
**CPU count:**  
* The number of cpus you would like to run per job at maximum in multithreaded steps (Mapping reads and Sorting BAMs).  
**Sample Filter File:**  
* The filter file is a space delimited text file with the name of the filter in question, the minimum value, the maximum value, and the step size:  

```
DPLT 1.0 3.0 0.2  
DPGT 0.2 1.0 0.2  
ABGT 0 0.5 0.05 
ABLT 0.5 0.95 0.05  
AD 0 5 1
GQ 5 95 5  
QUAL 5 200 5  
FS 0 30 5
QD 0 12 0.5
SOR 0 4 0.5
MQRankSumLT 0 2.5 0.5
MQRankSumGT -0.5 1.0 0.5
ReadPosRankSumLT 0 3 0.5
ReadPosRankSumGT -2.5 1.5 0.5
```

Available Filters are shown below. Idealized values are also put below for what a “gold standard” statistical value would be (ex. AB values should be around 0.5 for heterozygotes).  
For filters that have a GT or LT this refers to the filter either being “greater than” or “less than”(ex. DPGT establishes the lower bound recommendation for depth filters, ABLT establishes the upper bound recommendation for allele balance filters, etc.):  
**DPGT:** The lower bound of the depth filter, in scaled coverage (>=)  
**DPLT:** The upper bound of depth filter, in scaled coverage (<=)  
**ABGT:** Upper bound of proportion of ALT alleles to depth in child, 0.5 is ideal. (>=)  
**ABLT:** Lower bound of proportion of ALT alleles to depth in child, 0.5 is ideal. (<=)  
**AD:** Allele Depth of ALT allele in parents, 0 is ideal. (<=)  
**GQ:** Scaled Likelihood of genotype (>=)  
**QUAL:** Scaled likelihood of variation (>=)  
**QD:** Scaled Quality by Depth, Part of GATK Best Practices Hard Filter(>=)  
**FS:** Fisher Strand Test, checking for forward/reverse strand bias in heterozygotes. Part of GATK Best Practices Hard Filter, 0 is ideal. (<=)  
**SOR:** Strand Odds Ratio, checking for forward/reverse strange bias. part of GATK Best Practices Hard Filter. 0 is ideal. (<=)  
**MQRankSumGT:** Mapping Quality Rank Sum Test lower bound, compares mapping quality between the REF/ALT alleles. 0 is ideal. (>=)  
**MQRankSumLT:** Mapping Quality Rank Sum Test upper bound, compares mapping quality between the REF/ALT alleles. 0 is ideal. (<=)  
**ReadPosRankSumGT:** Read Position Rank Sum Test lower bound, compares whether the position of the variant on the REF/ALT reads are the same. 0 is ideal. (>=)  
**ReadPosRankSumLT:** Read Position Rank Sum Test upper bound, compares whether the position of the variant on the REF/ALT reads are the same. 0 is ideal. (<=)  

For a deeper explanation of each of the GATK Best Practices Hard Filter statistics, check here: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
###Demo Command/Job Submission:  
First, you must create the config file from which the workflow will read the user specifications:  
```
python DETECT/run_pipeline.py \
-R DETECT/demo/reference.fa \
-F DETECT/demo/filter_file.txt \
-U 2e-6 \
-O DETECT/demo/demo_workdir/best_filters.txt \
-V DETECT/demo/demo_variants.vcf \
--trio \
-P "dad,mom,junior" \
-C "10,20,30" \
-RL 100 -FL 300 -SD 30 \
-CL DETECT/demo/chrom_list.txt \
--cpus 12 -SP ~/app_list.txt -WD DETECT/demo/demo_workdir/  
```
Then, you can submit the snakemake job that will submit all subjobs. Note that this is more of a template, and the command may need to be altered to run on your cluster based on its SLURM configuration:   
```
sbatch -n1 --job-name demo_detect_superjob \
-o DETECT/demo/demo_workdir/demo_detect.out \
-e DETECT/demo/demo_workdir/demo_detect.err \
--wrap "snakemake -p --configfile  DETECT/demo/demo_workdir/config/config.json \
-s DETECT/Snakefile --default-resources mem_mb=8000 --scheduler greedy -j 100 \
--latency-wait 60 --keep-target-files --rerun-incomplete --cluster \
\" sbatch -n {threads} --mem={resources.mem_mb} -t 01:00:00 \
-o DETECT/demo/demo_workdir/logs/{rulename}.{jobid}.out \
-e DETECT/demo/demo_workdir/logs/{rulename}.{jobid}.err\" --forceall"
```
###My Job has run out of walltime!  
In the case that your DETECT job has run out of walltime, do not worry! Snakemake will pick up where it left off.  
Simply run this unlock command:
`snakemake --configfile <working_directory>/config/config.json -s <DETECT_directory>/DETECT/Snakefile --unlock`  

And then resubmit your job. It should continue from the last completed step. If an error occurred in a step, please put a support ticket into the repository, and I will be happy to help ASAP.  
###Best Practices:  
Once you have the filter recommendations, you may notice that the best filter is on one of the bounds of your parameter space (e.g. in the filter file above, you have a recommendation of 300 for QUAL). If this is the case, it is possible that there is a better filter beyond the bounds of your specifications, and thus it is recommended you expand your bounds and try again. Due to DETECT being implemented in snakemake, the workflow will only run the new filter thresholds.  

###Advanced Usage:
If you are familiar with the structure of JSON files, DETECT takes a JSON file as its input. For an example of an input file, see DETECT/demo/demo_workdit/config/config.json for a template guide. 
