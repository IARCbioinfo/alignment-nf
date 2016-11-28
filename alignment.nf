#! /usr/bin/env nextflow
// usage : ./alignment.nf --input_folder input/ --cpu 8 --mem 32 --ref hg19.fasta --RG "PL:ILLUMINA"
/*
vim: syntax=groovy
-*- mode: groovy;-*- */

// requirement:
// - samtools
// - samblaster
// - sambamba

//default values
params.help         = null
params.input_folder = '.'
params.fasta_ref    = 'hg19.fasta'
params.cpu          = 8
params.mem          = 32
params.RG           = ""
params.fastq_ext    = "fastq.gz"
params.suffix1      = "_1"
params.suffix2      = "_2"
params.out_folder   = "results_alignment"

if (params.help) {
    log.info ''
    log.info '-------------------------------------------------------------'
    log.info 'NEXTFLOW WHOLE EXOME/GENOME ALIGNMENT OR REALIGNMENT PIPELINE'
    log.info '-------------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run alignment.nf --input_folder input/ --fasta_ref hg19.fasta [--cpu 8] [--mem 32] [--RG "PL:ILLUMINA"] [--fastq_ext fastq.gz] [--suffix1 _1] [--suffix2 _2] [--out_folder output/]'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --input_folder   FOLDER                  Folder containing BAM or fastq files to be called.'
    log.info '    --fasta_ref          FILE                    Reference fasta file (with index).'
    log.info 'Optional arguments:'
    log.info '    --cpu          INTEGER                 Number of cpu used by bwa mem and sambamba (default: 8).'
    log.info '    --mem          INTEGER                 Size of memory used by sambamba (in GB) (default: 32).'
    log.info '    --RG           STRING                  Samtools read group specification with "\t" between fields.'
    log.info '                                           e.g. --RG "PL:ILLUMINA\tDS:custom_read_group".'
    log.info '                                           Default: "ID:bam_file_name\tSM:bam_file_name".'
    log.info '    --fastq_ext      STRING                Extension of fastq files (default : fastq.gz)'
    log.info '    --suffix1        STRING                Suffix of fastq files 1 (default : _1)'
    log.info '    --suffix2        STRING                Suffix of fastq files 2 (default : _2)'
    log.info '    --out_folder     STRING                Output folder (default: results_realignment).'
    log.info ''
    exit 1
}

//read files
fasta_ref     = file( params.fasta_ref )
fasta_ref_fai = file( params.fasta_ref+'.fai' )
fasta_ref_sa  = file( params.fasta_ref+'.sa' )
fasta_ref_bwt = file( params.fasta_ref+'.bwt' )
fasta_ref_ann = file( params.fasta_ref+'.ann' )
fasta_ref_amb = file( params.fasta_ref+'.amb' )
fasta_ref_pac = file( params.fasta_ref+'.pac' )

mode = 'fastq'
if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.fastq_ext}/ }.size() > 0){
    println "fastq files found, proceed with alignment"; files = Channel.fromPath( params.input_folder+'/*.${params.fastq_ext}' )
}else{
    if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0){
        println "BAM files found, proceed with realignment"; mode ='bam'; files = Channel.fromPath( params.input_folder+'/*.bam' )
    }else{
        println "ERROR: input folder contains no fastq nor BAM files"; System.exit(0)
    }
}

if(mode=='bam'){
    process bam_realignment {

        cpus params.cpu
        memory params.mem+'G'
        perJobMemLimit = true
  
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
        file('*realigned.bam*') into outputs

        shell:
        bam_tag = bam.baseName
        '''
        set -o pipefail
        samtools collate -uOn 128 !{bam_tag}.bam tmp_!{bam_tag} | samtools fastq - | bwa mem -M -t!{task.cpus} -R "@RG\\tID:!{bam_tag}\\tSM:!{bam_tag}\\t!{params.RG}" -p !{fasta_ref} - | samblaster --addMateTags | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t !{task.cpus} -m !{params.mem}G --tmpdir=!{bam_tag}_tmp -o !{bam_tag}_realigned.bam /dev/stdin
        '''
    }
}else if(mode=='fastq'){
    keys1 = file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.suffix1}${params.fastq_ext}/ }.collect { it.getName() }
                                                                                                               .collect { it.replace("${params.suffix1}${params.fastq_ext}",'') }
    keys2 = file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.suffix2}${params.fastq_ext}/ }.collect { it.getName() }
                                                                                                               .collect { it.replace("${params.suffix2}${params.fastq_ext}",'') }
    if ( !(keys1.containsAll(keys2)) || !(keys2.containsAll(keys1)) ) {println "\n ERROR : There is at least one fastq without its mate, please check your fastq files."; System.exit(0)}

    // Gather files ending with _1 suffix
    reads1 = Channel
    .fromPath( params.input_folder+'/*'+params.suffix1+params.fastq_ext )
    .map {  path -> [ path.name.replace("${params.suffix1}${params.fastq_ext}",""), path ] }

    // Gather files ending with _2 suffix
    reads2 = Channel
    .fromPath( params.input_folder+'/*'+params.suffix2+params.fastq_ext )
    .map {  path -> [ path.name.replace("${params.suffix2}${params.fastq_ext}",""), path ] }

    // Match the pairs on two channels having the same 'key' (name) and emit a new pair containing the expected files
    readPairs = reads1
    .phase(reads2)
    .map { pair1, pair2 -> [ pair1[1], pair2[1] ] }

    process fastq_alignment {

        cpus params.cpu
        memory params.mem+'GB'    
        tag { fastq_tag }
        publishDir params.out_folder, mode: 'move'

        input:
        file pair from readPairs
        file fasta_ref
        file fasta_ref_fai
        file fasta_ref_sa
        file fasta_ref_bwt
        file fasta_ref_ann
        file fasta_ref_amb
        file fasta_ref_pac
            
        output:
        //file('*aligned.bam*') into outputs

        shell:
        fastq_tag = fastq.baseName
        '''
        set -o pipefail
        bwa mem -M -t!{task.cpus} -R "@RG\\tID:!{fastq_tag}\\tSM:!{fastq_tag}\\t!{params.RG}" !{fasta_ref} !{pair[0]} !{pair[1]} | samblaster --addMateTags | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t !{task.cpus} -m !{params.mem}G --tmpdir=!{fastq_tag}_tmp -o !{fastq_tag}_aligned.bam /dev/stdin
        '''
    }
}
