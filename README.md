# DETECT (DNM Extraction Through Empirical Cutoff Thresholds)
DETECT is a simulation-based workflow that recommends filter thresholds that can be used to directly estimate mutation rate. By populating DNMs in a simulated trio at a specified mutation rate, we can determine filter thresholds that isolate as many DNMs as possible in the simulation, while also limiting the number of False Positives(FPs). These filter thresholds can then be applied to a real dataset with high confidence that you are isolating real mutations as well. DETECT has many options to tailor the simulation to real datasets including sequencing technology, coverage and variant datasets as input. DETECT has only been tested and is only functional on diploid, sexually reproducing organisms.

## Setting Up
### Environment Installation
```
git clone https://github.com/PfeiferLab/DETECT.git
conda env create -n DETECT -f DETECT/DETECT.yml  
source activate DETECT
```

## Quickstart
### Required Inputs
**Reference Genome:** Reference genome to be used in your real data workflow, where the simulated read data will come from. Must have a dictionary file (GATK CreateSequenceDictionary) and be bwa indexed (bwa index) e. g. reference.fa 

**Mutation Rate:** Number of mutations to populate the VCF(if >1), mutation rate to use to populate mutations(if <1), or positions file where mutations will be placed (positions file/VCF).

**Input Variants:** VCF file of trio containing variants to be used as False Positives. Must be without mendelian violations (MVs), and must be indexed (e.g. GATK IndexFeatureFile). 

**Read Length:** Length of the reads used in the real dataset. e. g. 100  

**Coverage:** Comma-delimited string of the coverages of the sire,dam,offspring. e. g. “30,40,50”  

**Output File:** Output file name of consolidated filter recommendations. 

**Pedigree:** Comma delimited string of the names of sire, dam, and offspring in the VCF. Required if trio VCF provided. (ex. “dad,mom,junior”)

### Optional Inputs

**Known Variants:** VCF file containing population-level variation to be used during BQSR. Default: None.

**Fragment Length:** Mean length of the fragment size distribution of the real data. Default: 300. (Mason Default)  

**Fragment Length Standard Deviation:** Standard deviation of the fragment size distribution of the real data. Default: 30 (Mason Default) 

**Chromosome list:** File of chromosome names to be simulated, one per line. Default: All contigs.  

**CPU count:** The number of cpus you would like to run per job at maximum in multithreaded steps (Mapping reads and Sorting BAMs).

**Number of Iterations:** The number of simulations you wish to perform. Default: 1

## Demo Command/Job Submission:
First, you must create the config file from which the workflow will read the user specifications:  
```
python DETECT/create_config.py \
-R DETECT/demo/reference.fa \
-U 2e-6 \
-O DETECT/demo/demo_workdir/best_filters.txt \
-V DETECT/demo/demo_variants.vcf \
-P "dad,mom,junior" \
-C "10,20,30" \
-RL 100 -FL 300 -SD 30 \
-CL DETECT/demo/chrom_list.txt \
--cpus 12 -WD DETECT/demo/demo_workdir/
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

## Understanding the Output File

When DETECT finishes, it creates one, tab delimited output file per iteration. The columns are:

**Filter**: Name of the filter.

**min/max**: denoting whether this is an upper or lower bound filter.
**Average:** Average value of summary statistic across de novo mutations(DNMs) that survived the pipeline.
**original_mutations:** the original number of DNMs populated into the simulation.
**total_sites:** the total number of sites in the final Mendelian Violation(MV) VCF.
**total_mutations:** the total number of DNMs that survived the pipeline.
**total_mutation_mut_recall:** recall of DNMs relative to the number of DNMs that were originally populated into the simulation.
**total_mutation_precision:** precision of DNMs relative to the total number of sites in the MV VCF.
**total_polymorphisms:** total number of sites that are present in the MV VCF that are miscalled polymorphisms.
**total_other_sites:** total number of sites that are present in the MV VCF that are neither DNMs nor miscalled polymorphisms.
**recommendation**: recommended filter values based on percentile cutoffs.
**filter_mutations:** number of DNMs retained after the recommended filter is applied to the MV VCF.
**filter_mutation_recall:** recall of DNMs after the recommended filter is applied to the MV VCF.
**filter_mutation_precision:** precision of DNMs after the recommended filter is applied to the MV VCF.
**filter_polymorphisms:** number of sites that are polymorphisms after the recommended filter is aplied to the MV VCF.
**filter_other_sites:** number of sites that are neither DNMs nor miscalled polymorphisms after the recommended filter is applied to the MV VCF.

The filternames are explained here:

| Filtername        | Full Name                  | Meaning                                                                           | Ideal Value                      | Types of filter | Filter applications   |
| ----------------- | -------------------------- | --------------------------------------------------------------------------------- | -------------------------------- | --------------- | --------------------- |
| DP                | Depth                      | Sequencing depth of site                                                          | average genomic depth            | min/max         | parents and offspring |
| GQ                | Genotype Quality           | Phred-scaled Genotype Quality of site                                             | as high as possible, maxed at 99 | min             | parents and offspring |
| QUAL              | Quality                    | QUAL score for presence of variation at the site                                  | as high as possible              | min             | per site              |
| AD                | Allele Depth               | number of reads with alternate alleles                                            | 0                                | max             | parents               |
| AB                | Allele Balance             | ratio of alternate reads to total depth of site                                   | 0.5                              | min/max         | offspring             |
| QD                | QualDepth                  | QUAL of site normalized by DP                                                     | as high as possible              | min             | per site              |
| MQRankSum         | Mapping Quality Rank Sum   | quantifies bias in mapping quality of reads that map to variant                   | 0                                | min/max         | per site              |
| ReadPosRankSum    | Read Position Rank Sum     | quantifies bias in the position of the variant on reads                           | 0                                | min/max         | per site              |
| FS                | Fisher Strand              | quantifies forward/reverse strand bias                                            | 0                                | max             | per site              |
| SOR               | Strand Odds Ratio          | quntifies forward/reverse strand bias                                             | 0                                | max             | per site              |
| parent.reassembly | parental reassembly filter | quantifies presence of reassembly during variant calling                          | 0                                | max             | parents               |
| child.reassembly  | child reassembly filter    | quantifies the difference in depth between the pre-calling BAM and reassembly BAM | 0                                | max             | offspring             |

For a deeper explanation of each of the GATK Best Practices Hard Filter statistics(QD and below on the table above), check here: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants  

For a closer look at the results per run, there are also output files in the run_outputs/ directory within your working directory. 
## My Job has run out of walltime!
In the case that your DETECT job has run out of walltime, do not worry! Snakemake will pick up where it left off.  
Simply run this unlock command:  
`snakemake --configfile <working_directory>/config/config.json -s <DETECT_directory>/DETECT/Snakefile --unlock`  

And then resubmit your job. It should continue from the last completed step. If an error occurred in a step, please put a support ticket into the repository, and I will be happy to help ASAP.   


## Advanced Usage:
If you are familiar with the structure of JSON files, DETECT takes a JSON file as its input. For an example of an input file, see DETECT/demo/demo_workdir/config/config.json for a template guide. 
