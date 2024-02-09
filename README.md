# alignment-nf

## Nextflow pipeline for BAM realignment or fastq alignment
[![CircleCI](https://circleci.com/gh/IARCbioinfo/alignment-nf.svg?style=svg)](https://circleci.com/gh/IARCbioinfo/alignment-nf)
[![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/r/iarcbioinfo/alignment-nf/)
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/4522)

![Workflow representation](WESpipeline.png?raw=true "Scheme of alignment/realignment Workflow")

## Description

Nextflow pipeline using BWA to perform BAM realignment or fastq alignment and QC, with/without base quality score recalibration.

## Dependencies

**A conda receipe, and docker and singularity containers are available with all the tools needed to run the pipeline (see "Usage")**

1. Nextflow : for common installation procedures see the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository.

#### Basic fastq alignment
2. [*bwa2*](https://github.com/bwa-mem2/bwa-mem2) (default) or [*bwa*](https://github.com/lh3/bwa)
3. [*samblaster*](https://github.com/GregoryFaust/samblaster)
4. [*sambamba*](https://github.com/lomereiter/sambamba)

#### BAM files realignment
5. [*samtools*](http://samtools.sourceforge.net/)

#### Adapter sequence trimming
6. [*AdapterRemoval*](https://github.com/MikkelSchubert/adapterremoval)

#### ALT contigs handling
7. The *k8* javascript execution shell (e.g., available in the [*bwakit*](https://sourceforge.net/projects/bio-bwa/files/bwakit/) archive); must be in the PATH
8. javascript bwa-postalt.js and the additional fasta reference *.alt* file from [*bwakit*](https://github.com/lh3/bwa/tree/master/bwakit) must be in the same directory as the reference genome file.

#### QC
9. [Qualimap](http://qualimap.bioinfo.cipf.es). 
10. [Multiqc](http://multiqc.info). 

#### Base quality score recalibration
11. [GATK4](https://software.broadinstitute.org/gatk/guide/quickstart); wrapper 'gatk' must be in the path
12. [GATK bundle](https://software.broadinstitute.org/gatk/download/bundle) VCF files with lists of indels and SNVs (recommended: Mills gold standard indels VCFs, dbsnp VCF), and corresponding tabix indexes (.tbi)
 
## Reference

In alignment processes, the reference FASTA file serves as a static template against which sequencing reads are aligned. It represents the backbone or framework of the genome against which variations and differences in individual sequenced DNA fragments are identified. The reference FASTA file is often obtained from publicly available databases, such as the NCBI, UCSC or GATK. 

| Name | Example value | Description |
|------|---------------|-------------| 
|--ref | hg38.fasta    | genome reference  with its index files |

## Indexes

To access the reference in a more efficient way, BWA need indexes files (*.fai*, *.sa*, *.bwt*, *.ann*, *.amb*, *.pac* and *.dict*) and bwa-mem2 also require (*.bwt.8bit.32 and *.0123). Indexes for the reference genome must be generated and stored in the same directory as the corresponding fasta file.

1. #### BWA indexes
```bash
bwa-mem index reference.fasta
bwa-mem2 index reference.fasta
```
2. #### Reference FAI Index for samtools :
```bash
samtools faidx reference.fasta
```
3. #### Reference Dictionary Index for GATK :
```bash
gatk CreateSequenceDictionary -R reference.fasta
```

If you wish to utilize tools within a container to generate indexes, you can proceed as follows:
```bash
singularity exec -B $(pwd):/mnt/ shub://IARCbioinfo/alignment-nf:v1.3 bwa-mem index reference.fasta
```

## Input (using a file)
| Type        | Description   |
|-------------|---------------|
|--input_file | Input file (comma-separated) with 4 columns: SM (sample name), RG (read group ID), pair1 (first fastq of the pair), and pair2 (second fastq of the pair). |

**You must provide the relative or complete path to the fastq files, not just the basename**
**Samples with the same name but different RG will be merged into one unique bam file.**

## Input (using a folder path)
| Type          | Description   |
|---------------|---------------|
|--input_folder | a folder with fastq files or bam files |

**With --input_folder, also use fastq_ext and suffix options to identify your pairs**

| Name | Default value  | Description |
|------------|----------|-------------| 
|--fastq_ext | fastq.gz | extension of fastq files|
|--suffix1   | \_1      | suffix for second element of read files pair|
|--suffix2   | \_2      | suffix for second element of read files pair|
|--RG        | PL:ILLUMINA | sequencing information for aligned (for *bwa*)|

## Usage
To run the pipeline with a file *input.txt* describing samples and fastq, one can type:
```bash
nextflow run iarcbioinfo/alignment-nf -r DLS2 -profile singularity --input_file input.txt --ref hg38.fasta --output_folder align
```
To run the pipeline on a series of fastq or BAM files in folder *input* and a fasta reference file hg38.fasta, one can type:
```bash
nextflow run iarcbioinfo/alignment-nf -r DLS2 -profile singularity --input_folder input/ --ref hg38.fasta --output_folder align
```
To run the pipeline without singularity just remove "-profile singularity". Alternatively, one can run the pipeline using a docker container (-profile docker), apptainer (-profile apptainer) or the conda receipe containing all required dependencies (-profile conda).

## Optional Parameters

| Name  | Description |
|-------|-------------| 
|--help | print usage and optional parameters |
|--bwa_option_M  | Trigger the -M option in bwa and the corresponding compatibility option in samblaster (marks shorter split hits as secondary) |

| Name      | Default value | Description |
|-----------|---------------|-------------| 
|--output_folder  | .  | Output folder for aligned BAMs|
|--cpu            | 8  | number of CPUs |
|--mem            | 32 | memory|

## Enable adapter trimming
To use the adapter trimming step, you must add the ***--trim* option**, as well as satisfy the requirements above mentionned. For example:
```bash
nextflow run iarcbioinfo/alignment-nf -r DLS2 -profile singularity --input_file input.txt --ref hg38.fasta --output_folder align --trim
```
| Name  | Description |
|-----------|-------------| 
|--trim     | enable adapter sequence trimming|

| Name      | Default value | Description |
|-----------|---------------|-------------| 
|--adapterremoval_opt | null | Command line options for AdapterRemoval |

## Enable ALT mode
To use the alternative contigs handling mode, you must provide the **path to an ALT aware genome reference** (e.g., hg38) AND add the ***--alt* option**, as well as satisfy the above-mentionned requirements. For example:
```bash
nextflow run iarcbioinfo/alignment-nf -r DLS2 -profile singularity  --input_file input.txt --ref hg38.fasta --output_folder align --alt
```
**hg38 and ALT aware genome references come with an additional file with the .alt extension that should be placed with the reference fasta and other indexes**

| Name  | Description |
|-------|-------------| 
|--alt  | enable alternative contig handling (for reference genome hg38)|

**If you are not using the container, you must also provide path to the bwa-postalt.js script from the bwa-kit : --postaltjs /user/bin/bwa-0.7.15/bwakit/bwa-postalt.js**

| Name       | Default value   | Description |
|------------|-----------------|-------------| 
|--postaltjs | bwa-postalt.js" | path to postalignment javascript *bwa-postalt.js*. (Already accessible in the container)|

## Enable base quality score recalibration
To use the base quality score recalibration step, you must provide the **path to 2 GATK bundle VCF files** with lists of known snps and indels, respectively, AND add the ***--recalibration* option**, as well as satisfy the requirements above mentionned. For example:
```bash
nextflow run iarcbioinfo/alignment-nf -r DSL2 -profile singularity  --input_file input.txt --ref hg38.fasta --output_folder align --snp_vcf GATKbundle/dbsnp.vcf.gz --indel_vcf GATKbundle/Mills_1000G_indels.vcf.gz --recalibration
```
| Name  | Description |
|-----------|-------------| 
|--recalibration  | perform quality score recalibration (GATK)|

| Name       | Default value | Description |
|------------|---------------|-------------| 
|--cpu\_BQSR | 2  | number of CPUs for GATK base quality score recalibration |
|--mem\_BQSR | 10 | memory for GATK base quality score recalibration |
|--snp_vcf   | dbsnp.vcf   | path to SNP VCF from GATK bundle (default : dbsnp.vcf) |
|--indel_vcf | Mills_1000G_indels.vcf | path to indel VCF from GATK bundle (default : Mills_1000G_indels.vcf) |

## Use bwa-mem instead of bwa-mem2
To utilize bwa-mem, one may input the following command:
```bash
nextflow run iarcbioinfo/alignment-nf -r DLS2 -profile singularity  --input_folder input/ --ref hg38.fasta --output_folder output --bwa_mem "bwa mem"
```
| Name       | Default value | Description |
|------------|---------------|-------------| 
|--bwa_mem   | bwa-mem2 mem  | bwa-mem command; use "bwa mem" to switch to regular bwa-mem (both are in the docker and singularity containers) |

## Optional parameters for QCs

| Name       | Default value | Description |
|------------|---------------|-------------| 
|--feature_file   | null |    Path to feature file for qualimap |
|--multiqc_config | null | config yaml file for multiqc |

## Output 
  | Type      | Description     |
  |-----------|---------------|
  | BAM/    | folder with BAM and BAI files of alignments or realignments |
  | QC/BAM/multiqc_qualimap_flagstat_*report.html  | multiQC report for qualimap and samtools flagstat (duplicates) |
  | QC/BAM/multiqc_qualimap_flagstat_*report_data  | data used for the multiQC report |
  | QC/qualimap/file_BQSRecalibrated.stats.txt  | qualimap summary file |
  | QC/qualimap/file_BQSRecalibrated/  | qualimap files |
  | QC/BAM/BQSR/  | GATK base quality score recalibration outputs (tables and pdf comparing scores before/after recalibration)|
  
## Directed Acyclic Graph
[![DAG](dag.png)](http://htmlpreview.github.io/?https://github.com/IARCbioinfo/alignment-nf/blob/dev/dag.html)

## FAQ
### Why did Indel realignment disappear from version 1.0?
Indel realignment was removed following new GATK best practices for pre-processing.

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------| 
  | Nicolas Alcala*    | AlcalaN@fellows.iarc.fr    | Developer to contact for support |
  | Catherine Voegele    |     VoegeleC@iarc.fr | Tester |
  | Vincent Cahais | CahaisV@iarc.fr | Developer |
  | Alexis Robitaille | RobitailleA@students.iarc.fr | Tester |
  
