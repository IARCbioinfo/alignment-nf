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
params.mem_sambamba = 1
params.RG           = "PL:ILLUMINA"
params.fastq_ext    = "fastq.gz"
params.suffix1      = "_1"
params.suffix2      = "_2"
params.out_folder   = "."
params.intervals    = ""
params.GATK_bundle       = "bundle"
params.GATK_folder  = "."
params.indel_realignment = "false"
params.recalibration = "false"
params.js           = "k8"
params.postaltjs    = "bwa-postalt.js"
params.alt          = "false"
params.trim         = "false"

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
    log.info '    --input_folder   FOLDER                  Folder containing BAM or fastq files to be aligned.'
    log.info '    --fasta_ref          FILE                    Reference fasta file (with index).'
    log.info 'Optional arguments:'
    log.info '    --indel_realignment                    Performs local indel realignment (default: no).'
    log.info '    --cpu          INTEGER                 Number of cpu used by bwa mem and sambamba (default: 8).'
    log.info '    --mem          INTEGER                 Size of memory used by sambamba (in GB) (default: 32).'
    log.info '    --RG           STRING                  Samtools read group specification with "\t" between fields.'
    log.info '                                           e.g. --RG "PL:ILLUMINA\tDS:custom_read_group".'
    log.info '                                           Default: "ID:bam_file_name\tSM:bam_file_name".'
    log.info '    --fastq_ext      STRING                Extension of fastq files (default : fastq.gz)'
    log.info '    --suffix1        STRING                Suffix of fastq files 1 (default : _1)'
    log.info '    --suffix2        STRING                Suffix of fastq files 2 (default : _2)'
    log.info '    --out_folder     STRING                Output folder (default: results_alignment).'
    log.info ''
    exit 1
}

//read files
fasta_ref = file(params.fasta_ref)
fasta_ref_fai = file( params.fasta_ref+'.fai' )
fasta_ref_sa = file( params.fasta_ref+'.sa' )
fasta_ref_bwt = file( params.fasta_ref+'.bwt' )
fasta_ref_ann = file( params.fasta_ref+'.ann' )
fasta_ref_amb = file( params.fasta_ref+'.amb' )
fasta_ref_pac = file( params.fasta_ref+'.pac' )
fasta_ref_alt = file( params.fasta_ref+'.alt' )

mode = 'fastq'
if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.fastq_ext}/ }.size() > 0){
    println "fastq files found, proceed with alignment"
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
        tag { file_tag }
        
        input:
        file infile from files
	file fasta_ref
	file fasta_ref_fai
	file fasta_ref_sa
	file fasta_ref_bwt
	file fasta_ref_ann
	file fasta_ref_amb
	file fasta_ref_pac
	file fasta_ref_alt
     
        output:
	set val(file_tag_new), file("${file_tag_new}.bam") into bam_files
	file("${file_tag_new}.bam.bai") into bai_files
	if( (params.recalibration=="false")&(params.indel_realignment=="false") ) publishDir params.out_folder, mode: 'move'

        shell:
	file_tag = infile.baseName
	file_tag_new=file_tag+'_realigned'
	if(params.trim!="false") file_tag_new=file_tag_new+'_trimmed'
	if(params.alt!="false")  file_tag_new=file_tag_new+'_alt'
	
	if(params.alt=="false"){
	  ignorealt='-j'
	  postalt=''
	}else{
	  ignorealt=''
	  postalt=params.js+' '+params.postaltjs+' '+params.fasta_ref+'.alt |'
	}
	if(params.trim=="false"){
	  preproc=''
	}else{
	  	
	  preproc='AdapterRemoval --interleaved --file1 /dev/stdin --output1 /dev/stdout |'
	}
	
        '''
        set -o pipefail
        samtools collate -uOn 128 !{file_tag}.bam tmp_!{file_tag} | samtools fastq - | !{preproc} bwa mem !{ignorealt} -M -t!{task.cpus} -R "@RG\\tID:!{file_tag}\\tSM:!{file_tag}\\t!{params.RG}" -p !{params.fasta_ref} - | !{postalt} samblaster --addMateTags | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t !{task.cpus} -m !{params.mem_sambamba}G --tmpdir=!{file_tag}_tmp -o !{file_tag_new}.bam /dev/stdin
        '''
    }
}
if(mode=='fastq'){
    println "fastq mode"
    keys1 = file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.suffix1}.${params.fastq_ext}/ }.collect { it.getName() }
                                                                                                               .collect { it.replace("${params.suffix1}.${params.fastq_ext}",'') }
    keys2 = file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.suffix2}.${params.fastq_ext}/ }.collect { it.getName() }
                                                                                                               .collect { it.replace("${params.suffix2}.${params.fastq_ext}",'') }
    if ( !(keys1.containsAll(keys2)) || !(keys2.containsAll(keys1)) ) {println "\n ERROR : There is at least one fastq without its mate, please check your fastq files."; System.exit(0)}

    println keys1

    // Gather files ending with _1 suffix
    reads1 = Channel
    .fromPath( params.input_folder+'/*'+params.suffix1+'.'+params.fastq_ext )
    .map {  path -> [ path.name.replace("${params.suffix1}.${params.fastq_ext}",""), path ] }

    // Gather files ending with _2 suffix
    reads2 = Channel
    .fromPath( params.input_folder+'/*'+params.suffix2+'.'+params.fastq_ext )
    .map {  path -> [ path.name.replace("${params.suffix2}.${params.fastq_ext}",""), path ] }

    // Match the pairs on two channels having the same 'key' (name) and emit a new pair containing the expected files
    readPairs = reads1
    .phase(reads2)
    .map { pair1, pair2 -> [ pair1[1], pair2[1] ] }

    println reads1
        
    process fastq_alignment {

        cpus params.cpu
        memory params.mem+'GB'    
        tag { file_tag }
        
        input:
        file pair from readPairs
	file fasta_ref
	file fasta_ref_fai
	file fasta_ref_sa
	file fasta_ref_bwt
	file fasta_ref_ann
	file fasta_ref_amb
	file fasta_ref_pac
	file fasta_ref_alt
                 
        output:
	set val(file_tag_new), file("${file_tag_new}.bam") into bam_files
	file("${file_tag_new}.bam.bai") into bai_files
	if( (params.recalibration=="false")&(params.indel_realignment=="false") ) publishDir params.out_folder, mode: 'move'

        shell:
        file_tag = pair[0].name.replace("${params.suffix1}.${params.fastq_ext}","")
	file_tag_new=file_tag
	if(params.trim!="false") file_tag_new=file_tag_new+'_trimmed'
	if(params.alt!="false")  file_tag_new=file_tag_new+'_alt'
	if(params.alt=="false"){
	  ignorealt='-j'
	  postalt=''
	}else{
	  ignorealt=''
	  postalt=params.js+' '+params.postaltjs+' '+params.fasta_ref+'.alt |'
	}
	if(params.trim=="false"){
		'''
        	set -o pipefail
		bwa mem !{ignorealt} -M -t!{task.cpus} -R "@RG\\tID:!{file_tag}\\tSM:!{file_tag}\\t!{params.RG}" !{params.fasta_ref} !{pair[0]} !{pair[1]} | !{postalt} samblaster --addMateTags | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t !{task.cpus} -m !{params.mem}G --tmpdir=!{file_tag}_tmp -o !{file_tag_new}.bam /dev/stdin
		'''
 	}else{
		'''
        	set -o pipefail
		AdapterRemoval --file1 !{pair[0]} --file2 !{pair[1]} --interleaved-output --output1 /dev/stdout | bwa mem !{ignorealt} -M -t!{task.cpus} -R "@RG\\tID:!{file_tag}\\tSM:!{file_tag}\\t!{params.RG}" -p !{params.fasta_ref} - | !{postalt} samblaster --addMateTags | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t !{task.cpus} -m !{params.mem}G --tmpdir=!{file_tag}_tmp -o !{file_tag_new}.bam /dev/stdin
		'''
	}
	
     }
}

if(params.indel_realignment != "false"){
        // Local realignment around indels
        process indel_realignment {
            cpus params.cpu
            memory params.mem+'G'
            tag { file_tag }
            input:
	    set val(file_tag), file("${file_tag}.bam") from bam_files
	    file("${file_tag}.bam.bai") from bai_files
            output:
            file("${file_tag_new}_target_intervals.list") into indel_realign_target_files
            set val(file_tag_new), file("${file_tag_new}.bam") into bam_files2
	    file("${file_tag_new}.bam.bai") into bai_files2
	    if(params.recalibration=="false") publishDir params.out_folder, mode: 'move'
	    
            shell:
	    file_tag_new=file_tag+'_indelrealigned'
            '''
	    indelsvcf=`ls !{params.GATK_bundle}/*indels*.vcf* | grep -v ".tbi"`
	    knowncom=''
	    for ll in $indelsvcf; do knowncom=$knowncom' -known '$ll; done
            java -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt !{params.cpu} -R !{params.fasta_ref} -I !{file_tag}.bam $knowncom -o !{file_tag_new}_target_intervals.list
            java -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T IndelRealigner -R !{params.fasta_ref} -I !{file_tag}.bam -targetIntervals !{file_tag_new}_target_intervals.list $knowncom -o !{file_tag_new}.bam
	    mv !{file_tag_new}.bai !{file_tag_new}.bam.bai
            '''
        }
}else{
    if(params.recalibration!= "false"){
    process no_indel_realignment {
        cpus '1'
        memory '100M'
        tag { file_tag }
        
        input:
        set val(file_tag), file("${file_tag}.bam") from bam_files
	file("${file_tag}.bam.bai") from bai_files
        output:
        set val(file_tag), file("${file_tag}.bam") into bam_files2
	file("${file_tag}.bam.bai") into bai_files2
	shell:
	'''
        touch !{file_tag}.bam
	'''
    }
    }
}

if(params.recalibration!= "false"){
// base quality score recalibration
   process base_quality_score_recalibration {
    cpus params.cpu
    memory params.mem+'G'
    tag { file_tag }
        
    input:
    set val(file_tag), file("${file_tag}.bam") from bam_files2
    file("${file_tag}.bam.bai") from bai_files2
    output:
    file("${file_tag_new}_recal.table") into recal_table_files
    file("${file_tag_new}_post_recal.table") into recal_table_post_files
    file("${file_tag_new}_recalibration_plots.pdf") into recal_plots_files
    set val(file_tag_new), file("${file_tag_new}.bam") into recal_bam_files
    file("${file_tag_new}.bam.bai") into recal_bai_files
    publishDir params.out_folder, mode: 'move'

    shell:
    file_tag_new=file_tag+'_BQSrecalibrated'
    '''
    indelsvcf=(`ls !{params.GATK_bundle}/*indels*.vcf* | grep -v ".tbi"`)
    dbsnpvcfs=(`ls !{params.GATK_bundle}/*dbsnp*.vcf* | grep -v ".tbi"`)
    dbsnpvcf=${dbsnpvcfs[@]:(-1)}
    knownSitescom=''
    for ll in $indelsvcf; do knownSitescom=$knownSitescom' -knownSites '$ll; done
    knownSitescom=$knownSitescom' -knownSites '$dbsnpvcf
    java -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T BaseRecalibrator -nct !{params.cpu} -R !{params.fasta_ref} -I !{file_tag}.bam $knownSitescom -L !{params.intervals} -o !{file_tag_new}_recal.table
    java -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T BaseRecalibrator -nct !{params.cpu} -R !{params.fasta_ref} -I !{file_tag}.bam $knownSitescom -BQSR !{file_tag_new}_recal.table -L !{params.intervals} -o !{file_tag_new}_post_recal.table		
    java -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T AnalyzeCovariates -R !{params.fasta_ref} -before !{file_tag_new}_recal.table -after !{file_tag_new}_post_recal.table -plots !{file_tag_new}_recalibration_plots.pdf	
    java -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T PrintReads -nct !{params.cpu} -R !{params.fasta_ref} -I !{file_tag}.bam -BQSR !{file_tag_new}_recal.table -L !{params.intervals} -o !{file_tag_new}.bam
    mv !{file_tag_new}.bai !{file_tag_new}.bam.bai
    '''
    }
}
