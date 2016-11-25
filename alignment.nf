#! /usr/bin/env nextflow
// usage : ./bam_realignment.nf --input_folder BAM/ --cpu 8 --mem 32 --ref hg19.fasta --RG "PL:ILLUMINA"
/*
vim: syntax=groovy
-*- mode: groovy;-*- */

if (params.help) {
    log.info ''
    log.info '--------------------------------------------------'
    log.info 'NEXTFLOW ALIGNMENT OR REALIGNMENT PIPELINE'
    log.info '--------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run bam_realignment.nf --input_folder BAM/ --cpu 8 --mem 32 --fasta_ref hg19.fasta'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --input_folder   FOLDER                  Folder containing BAM or fastq files to be called.'
    log.info '    --ref          FILE                    Reference fasta file (with index).'
    log.info 'Optional arguments:'
    log.info '    --cpu          INTEGER                 Number of cpu used by bwa mem and sambamba (default: 8).'
    log.info '    --mem          INTEGER                 Size of memory used by sambamba (in GB) (default: 32).'
    log.info '    --RG           STRING                  Samtools read group specification with "\t" between fields.'
    log.info '                                           e.g. --RG "PL:ILLUMINA\tDS:custom_read_group".'
    log.info '                                           Default: "ID:bam_file_name\tSM:bam_file_name".'
    log.info '    --out_folder   STRING                  Output folder (default: results_realignment).'
    log.info ''
    exit 1
}

//default values
params.input_folder = '.'
params.ref          = 'hg19.fasta'
params.cpu          = 8
params.mem          = 32
params.RG           = ""
params.out_folder   = "results_alignment"

//read files
fasta_ref     = file( params.fasta_ref )
fasta_ref_fai = file( params.fasta_ref+'.fai' )
fasta_ref_sa  = file( params.fasta_ref+'.sa' )
fasta_ref_bwt = file( params.fasta_ref+'.bwt' )
fasta_ref_ann = file( params.fasta_ref+'.ann' )
fasta_ref_amb = file( params.fasta_ref+'.amb' )
fasta_ref_pac = file( params.fasta_ref+'.pac' )

mode  = 'fastq'
files = Channel.fromPath( params.input_folder+'/*.fastq' )
              .ifEmpty { mode='bam'}
if(mode=='bam'){
        files = Channel.fromPath( params.input_folder+'/*.bam' )
                       .ifEmpty { error "Cannot find any fastq or bam file in: ${params.input_folder}" }
        }
              
process bam_realignment {

    cpus params.cpu
    memory params.mem+'GB'    
  
    tag { bam_tag }

    publishDir params.out_folder, mode: 'move'

    input:
    file bam from files
    file fasta_ref
    file fasta_ref_fai
    file fasta_ref_sa
    file fasta_ref_bwt
    file fasta_ref_ann
    file fasta_ref_amb
    file fasta_ref_pac

    output:
    //file('*realigned.bam*') into outputs

    shell:
    bam_tag = bam.baseName
    echo '''
    set -o pipefail
    samtools collate -uOn 128 !{bam_tag}.bam tmp_!{bam_tag} | samtools fastq - | bwa mem -M -t!{task.cpus} -R "@RG\\tID:!{bam_tag}\\tSM:!{bam_tag}\\t!{params.RG}" -p !{fasta_ref} - | samblaster --addMateTags | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t !{task.cpus} -m !{params.mem}G --tmpdir=!{bam_tag}_tmp -o !{bam_tag}_realigned.bam /dev/stdin
    '''
}



