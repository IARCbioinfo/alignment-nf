# alignment-nf

## Nextflow pipeline for BAM realignment or fastq alignment

![Workflow representation](WESpipeline.png?raw=true "Scheme of alignment/realignment Workflow")

## Description

Nextflow pipeline to perform BAM realignment or fastq alignment, with/without local indel realignment and base quality score recalibration.

## Dependencies

1. Nextflow : for commom installation procedures see the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository.

### Basic fastq alignment
2. [*bwa*](https://github.com/lh3/bwa)
3. [*samblaster*](https://github.com/GregoryFaust/samblaster)
4. [*sambamba*](https://github.com/lomereiter/sambamba)

### BAM files realignment
5. [*samtools*](http://samtools.sourceforge.net/)

### Adapter sequence trimming
6. [*AdapterRemoval*](https://github.com/MikkelSchubert/adapterremoval)

### ALT contigs handling
7. the *k8* javascript execution shell (e.g., available in the [*bwakit*](https://sourceforge.net/projects/bio-bwa/files/bwakit/) archive)
8. javascript bwa-postalt.js and the additional fasta reference *.alt* file from [*bwakit*](https://github.com/lh3/bwa/tree/master/bwakit) must be in the same directory as the reference genome file.

### Indel realignment
9. GATK [*GenomeAnalysisTK.jar*](https://software.broadinstitute.org/gatk/guide/quickstart)
10. [GATK bundle](https://software.broadinstitute.org/gatk/download/bundle) VCF files with lists of indels (recommended: 1000 genomes and Mills gold standard VCFs)

### Base quality score recalibration
11. GATK [*GenomeAnalysisTK.jar*](https://software.broadinstitute.org/gatk/guide/quickstart)
12. [GATK bundle](https://software.broadinstitute.org/gatk/download/bundle) VCF files with lists of indels and SNVs (recommended: 1000 genomes indels, Mills gold standard indels VCFs, dbsnp VCF)

To avoid installing the previous tools, install Docker. Docker installation is described in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository.

## Input 
 | Type      | Description     |
  |-----------|---------------|
  | --input_folder    | a folder with fastq files or bam files |

## Parameters

* #### Mandatory

| Name | Example value | Description |
|-----------|--------------|-------------| 
|--ref    | hg19.fasta | genome reference  with its index files (*.fai*, *.sa*, *.bwt*, *.ann*, *.amb*, *.pac*; in the same directory) |
|--output_folder   | . | Output folder for aligned BAMs|

* #### Optional

| Name | Default value | Description |
|-----------|--------------|-------------| 
|--cpu          | 8 | number of CPUs |
|--mem         | 32 | memory|
|--mem\_sambamba | 1 | memory for software *sambamba*|
|--RG           | PL:ILLUMINA | sequencing information for aligned (for *bwa*)|
|--fastq_ext    | fastq.gz | extension of fastq files|
|--suffix1      | \_1 | suffix for second element of read files pair|
|--suffix2      | \_2 | suffix for second element of read files pair|
|--bed    | | bed file with interval list|
|--GATK_bundle  | bundle | path to GATK bundle files|
|--GATK_folder  | . | path to GATK *GenomeAnalysisTK.jar* file |
|--js           | k8 | path to javascript interpreter *k8*|
|--postaltjs    | bwa-postalt.js" | path to postalignment javascript *bwa-postalt.js*|


* #### Flags

Flags are special parameters without value.

| Name  | Description |
|-----------|-------------| 
| --help | print usage and optional parameters |
|--trim     | enable adapter sequence trimming|
|--indel_realignment  | perform local indel realignment (GATK)|
|--recalibration  | perform quality score recalibration (GATK)|
|--alt         | enable alternative contig handling (for reference genome hg38)|

## Usage
```bash
nextflow run iarcbioinfo/alignment-nf --input_folder input --fasta_ref hg19.fasta --out_folder output
```
### Enable adapter trimming
To use the adapter trimming step, you must add the ***--trim* option**, as well as satisfy the requirements above mentionned. For example:
```bash
nextflow run iarcbioinfo/alignment-nf --input_folder input --fasta_ref reference/hs38DH.fa -out_folder output --trim
```

### Enable ALT mode
To use the alternative contigs handling mode, you must provide the **path to an ALT aware genome reference** (e.g., hg38) AND add the ***--alt* option**, as well as satisfy the requirements above mentionned. For example:
```bash
nextflow run iarcbioinfo/alignment-nf --input_folder input --fasta_ref reference/hs38DH.fa --js /user/bin/k8/k8 --postaltjs /user/bin/bwa-0.7.15/bwakit/bwa-postalt.js -out_folder output --alt
```
### Enable local indel realignment
To use the local indel realignment step, you must provide the **path to the GATK jar file**, the **GATK bundle folder**, AND add the ***--indel_realignment* option**, as well as satisfy the requirements above mentionned. For example:
```bash
nextflow run iarcbioinfo/alignment-nf --GATK_bundle GATKbundle/hg19 --input_folder input --fasta_ref reference/hg19.fa --GATK_folder /user/bin7GATK-3.6-0 --out_folder output --indel_realignment
```

### Enable base quality score recalibration
To use the base quality score recalibration step, you must provide the **path to the GATK jar file**, the **GATK bundle folder**, a **bed file**, AND add the ***--recalibration* option**, as well as satisfy the requirements above mentionned. For example:
```bash
nextflow run iarcbioinfo/alignment-nf --GATK_bundle GATKbundle/hg19 --input_folder input --fasta_ref reference/hg19.fa --GATK_folder /user/bin7GATK-3.6-0 --intervals reference/hg19_intervals.bed --out_folder output --recalibration
```

## Output 
  | Type      | Description     |
  |-----------|---------------|
  | output1    | ...... |
  | output2    | ...... |

## Directed Acyclic Graph


## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------| 
  | Nicolas Alcala*    | AlcalaN@fellows.iarc.fr    | Developer to contact for support |
  | Catherine Voegele    |            xx | Developer |
  | contrib3    |            xx | Tester |
  