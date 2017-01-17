# alignment-nf

## Overview of pipeline workflow

![workflow](WESpipeline.png?raw=true "Scheme of alignment/realignment Workflow")

## Prerequisites
For the basic fastq alignment without indel realignment, base quality score recalibration, nor alternative contif handling, the script requires the following programs:
- *bwa*
- *samblaster*
- *sambamba*
and the following files:
- a fasta reference genome with its index files (*.fai*, *.sa*, *.bwt*, *.ann*, *.amb*, *.pac*; in the same directory)
- a folder with fastq files

In addition, for the bam file realignment:
- *samtools*

For the ALT contigs handling, additional softwares and scripts are required:
- the *k8* javascript execution shell
- javascript bwa-postalt.js from *bwakit*
and also the additional *.alt* file must be in the same directory as the reference genome file.

For the indel realignment:
- GATK *GenomeAnalysisTK.jar*
- GATK bundle VCF files with lists of indels (recommended: 1000 genomes and Mills gold standard VCFs)

For the base quality score recalibration:
- GATK *GenomeAnalysisTK.jar*
- GATK bundle VCF files with lists of indels and SNVs (recommended: 1000 genomes indels, Mills gold standard indels VCFs, dbsnp VCF)
- bed file with intervals to be considered

## Usage
```bash
nextflow run alignment.nf --input_folder input --fasta_ref hg19.fasta --out_folder output
```
### Enable ALT mode
To use the alternative contigs handling mode, you must provide the path to an ALT aware genome reference (e.g., hg38) AND add the *--alt* option, as well as satisfy the requirements above mentionned. For example:
```bash
nextflow run alignment-nf --input_folder input --fasta_ref reference/hs38DH.fa --js /user/bin/k8/k8 --postaltjs /user/bin/bwa-0.7.15/bwakit/bwa-postalt.js -out_folder output --alt
```
### Enable local indel realignment
To use the local indel realignment step, you must provide the path to the GATK jar file, the GATK bundle folder AND add the *--indel_realignment* option, as well as satisfy the requirements above mentionned. For example:
```bash
nextflow run alignment-nf --GATK_bundle GATKbundle/hg19 --input_folder input --fasta_ref reference/hg19.fa --GATK_folder /user/bin7GATK-3.6-0 --out_folder output --indel_realignment
```

### Enable base quality score recalibration
To use the base quality score recalibration step, you must provide the path to the GATK jar file, the GATK bundle folder AND add the *--indel_realignment* option, as well as satisfy the requirements above mentionned. For example:
```bash
nextflow run alignment-nf --GATK_bundle GATKbundle/hg19 --input_folder input --fasta_ref reference/hg19.fa --GATK_folder /user/bin7GATK-3.6-0 --intervals reference/hg19_intervals.bed --out_folder output --recalibration
```

## All parameters
| **PARAMETER** | **DEFAULT** | **DESCRIPTION** |
|-----------|--------------:|-------------| 
| *--help* | null | print usage and optional parameters |
*--input_folder* | . | input folder |
*--fasta_ref*    | hg19.fasta | genome reference |
*--cpu*          | 8 | number of CPUs |
*--mem*         | 32 | memory|
*--mem\_sambamba* | 1 | memory for software *sambamba*|
*--RG*           | PL:ILLUMINA | sequencing information for aligned (for *bwa*)|
*--fastq_ext*    | fastq.gz | extension of fastq files|
*--suffix1*      | \_1 | suffix for second element of read files pair|
*--suffix2*      | \_2 | suffix for second element of read files pair|
*--out_folder*   | . | output folder for aligned BAMs|
*--intervals*    | | bed file with interval list|
*--GATK_bundle*  | bundle | path to GATK bundle files|
*--GATK_folder*  | . | path to GATK *GenomeAnalysisTK.jar* file |
*--indel_realignment* | false | perform local indel realignment (GATK)|
*--recalibration* | false | perform quality score recalibration (GATK)|
*--js*           | k8 | path to javascript interpreter *k8*|
*--postaltjs*    | bwa-postalt.js" | path to postalignment javascript *bwa-postalt.js*|
*--alt*          | false | enable alternative contig handling (for reference genome hg38)|
