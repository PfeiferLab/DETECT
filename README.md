## DETECT: a simulation framework optimizing <i>de novo</i> mutation detection across species and study designs

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


### Advanced Usage

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


## Understanding DETECT's Output

By default, DETECT evaluates a variety of computational filter criteria and thresholds to obtain species- and study design-specific best practice recommendations.

| filter criterion* | filter description                                                                               | bounds  | applied to           |
| ----------------- | ------------------------------------------------------------------------------------------------ | ------- | -------------------- |
| ```AD```                | allele depth                                                                               | max     | sire, dam            |
| ```AB```                | allele balance                                                                             | min/max | offspring            |
| ```DP```                | depth of coverage                                                                          | min/max | sire, dam, offspring |
| ```GQ```                | genotype quality                                                                           | min     | sire, dam, offspring |
| ```QUAL```              | quality score                                                                              | min     | sire, dam, offspring |
| ```parent.reassembly``` | scaled depth of coverage in the reassembled region during variant calling in the parent    | max     | sire, dam            |
| ```child.reassembly```  | scaled depth of coverage in the reassembled region during variant calling in the offspring | max     | offspring            |

<br>

DETECT creates one tab-delimited output file per iteration containing the recommended filter criteria and thresholds. The columns are:

* **filter**: Filter criterion.

* **min/max**: Denoting whether the recommended filter is an upper or lower bound filter.

* **average:** Average value of across identified de novo mutations (DNMs).

* **original_mutations:** The number of DNMs introduced in the offspring during the simulation.

* **total_sites:** The number of sites with Mendelian violations in the patterns of inheritance.

* **total_mutations:** The number of DNMs passing the recommended filter criterion and threshold.

* **total_mutation_mut_recall:** Recall (the number DNMs passing the recommended filter criterion and threshold divided by the number of DNMs introduced in the offspring during the simulation).

* **total_mutation_precision:** Precision (the number DNMs passing the recommended filter criterion and threshold divided by the number of sites with Mendelian violations in the patterns of inheritance).

* **total_polymorphisms:** The number of segregating variants with Mendelian violations in the patterns of inheritance (false positives).

* **total_other_sites:** The number of miscalled DNMs resulting from technical artefacts (false positives).

* **recommendation**: Recommended filter threshold based on the percentile cutoff.

* **filter_mutations:** The number of DNMs passing all recommended filter criteria and thresholds.

* **filter_mutation_recall:** Recall (the number DNMs passing all recommended filter criteria and thresholds divided by the number of DNMs introduced in the offspring during the simulation).

* **filter_mutation_precision:** Precision (the number DNMs passing all recommended filter criteria and thresholds divided by the number of sites with Mendelian violations in the patterns of inheritance).

* **filter_polymorphisms:** The number of segregating variants with Mendelian violations in the patterns of inheritance remaining after the application of all recommended filter criteria and thresholds (false positives).

* **filter_other_sites:** The number of miscalled DNMs resulting from technical artefacts remaining after the application of all recommended filter criteria and thresholds (false positives).

<i>For a closer look at the results per run, please refer to the output files in the run_outputs/ directory.</i> 

<br>

## My Job has run out of walltime!
In the case that your DETECT job has run out of walltime, do not worry! Snakemake will pick up where it left off.  
Simply run this unlock command:  
`snakemake --configfile <working_directory>/config/config.json -s <DETECT_directory>/DETECT/Snakefile --unlock`  
and then resubmit your job. It should continue from the last completed step.
