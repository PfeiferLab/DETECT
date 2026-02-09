# DETECT (DNM Extraction Through Empirical Cutoff Thresholds)

![DETECT logo](https://github.com/user-attachments/assets/5b1d31b1-84eb-4436-aa3b-c3d38e8e5be6)

DETECT allows researchers to obtained best practice recommendations for computational filter criteria and thresholds mitigating artefacts in <i>de novo</i> mutation (DNM) detection tailored to their specific study design and species of interest.


## Setting Up
### Environment Installation
```
git clone https://github.com/PfeiferLab/DETECT.git
mamba env create -n DETECT -f DETECT/DETECT.yml  
source activate DETECT
```

## Quickstart
### Inputs (required)

* **Reference Genome:** A reference assembly (.fasta) from which paired-end reads will be simulated. Must have an associated sequence dictionary file (.dict) and BWA index files (.amb, .ann, .bwt, .pac, .sa). 

* **Segregating Variants:** Segregating variants (.vcf). Must be without Mendelian violations and must be indexed (.idx).

* **Pedigree:** Comma delimited string of sample identifiers and genetic relationships of each individual in the trio (e.g., “sire,dam,offspring“). Required if trio data (.vcf) is provided. If no filial information is provided, DETECT will “generate” an offspring based on either the parental haplotypes (if parental information is available), or two user-specified, or randomly selected, haplotypes sampled from the population (if parental information is unavailable).

* **Coverage:** Comma-delimited string of the depth of coverages for the sire, dam, and offspring to be simulated (e.g., “30,40,50”)  

* **Read Length:** Length of the reads to be simulated (e.g., “150“)

* **Mutation Rate:** Either a rate (or number) of DNMs that should be introduced at random into the simulated reads of the offspring, or a file (in .bed or .vcf format) containing the positions where DNMs should be spiked in.


### Inputs (optional)

* **Chromosome list:** A list of chromosomes from which paired-end reads will be simulated (default: all).  

* **Known Variants:** Experimentally validated variants to facilitate base quality score recalibration (BQSR) of the simulated reads (default: none).

* **Fragment Length:** Mean length of the fragment size distribution (default: 300).  

* **Fragment Length Standard Deviation:** Standard deviation of the fragment size distribution (default: 30). 

* **Number of Iterations:** The number of replicates performed within a single run (default: 1).

* **Number of CPUs:** The number of CPUs per job (default: 1).


### Output

* **Output File:** Best practice recommendations for computational filter criteria and thresholds.


## Advanced Usage

* If you are familiar with the structure of JSON files, DETECT takes a JSON file as its input. For an example of an input file, see DETECT/demo/demo_workdir/config/config.json for a template guide. 
<br>


## Demo Command/Job Submission:
First, create a config file from which the workflow will read the user specifications:  
```
python DETECT/create_config.py \
-R DETECT/demo/reference.fa \
-U 2e-6 \
-O DETECT/demo/demo_workdir \
-V DETECT/demo/demo_variants.vcf \
-P "dad,mom,junior" \
-C "10,20,30" \
-RL 100 -FL 300 -SD 30 \
-CL DETECT/demo/chrom_list.txt \
--cpus 12 \
-WD DETECT/demo/demo_workdir/
```

Then, submit the Snakemake job that will submit all sub-jobs (note that the command may need to be altered according to the specific cluster environment / configuration):   
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
<br>


## Understanding the Output File

When DETECT finishes, it creates one, tab delimited output file per iteration. The columns are:

**Filter**: Name of the filter.

**min/max**: denoting whether this is an upper or lower bound filter.

**average:** Average value of summary statistic across de novo mutations(DNMs) that survived the pipeline.

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

<br>

## My Job has run out of walltime!
In the case that your DETECT job has run out of walltime, do not worry! Snakemake will pick up where it left off.  
Simply run this unlock command:  
`snakemake --configfile <working_directory>/config/config.json -s <DETECT_directory>/DETECT/Snakefile --unlock`  
and then resubmit your job. It should continue from the last completed step.
