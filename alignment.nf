#! /usr/bin/env nextflow

/*vim: syntax=groovy -*- mode: groovy;-*- */

// Copyright (C) 2017 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.input_folder = '.'
params.input_file   = null
params.ref          = 'hg19.fasta'
params.cpu          = 8
params.mem          = 32
params.RG           = "PL:ILLUMINA"
params.fastq_ext    = "fastq.gz"
params.output_folder = "."
params.bed          = ""
params.GATK_bundle  = "bundle"
params.GATK_folder  = "."
params.postaltjs    = "bwa-postalt.js"
params.indel_realignment = null
params.recalibration = null
params.help         = null
params.alt          = null
params.trim         = null


log.info ""
log.info "--------------------------------------------------------"
log.info "  alignment-nf 2.0.0: alignment/realignment workflow for whole exome/whole genome sequencing "
log.info "--------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""


if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ''
    log.info 'nextflow run iarcbioinfo/alignment.nf [-with-docker] --input_folder input/ --ref hg19.fasta [OPTIONS]'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '--input_folder   FOLDER                  Folder containing BAM or fastq files to be aligned.'
    log.info '--ref          FILE                    Reference fasta file (with index).'
    log.info '--output_folder     STRING                Output folder (default: results_alignment).'
    log.info ""
    log.info 'Optional arguments:'
    log.info '    --cpu          INTEGER                 Number of cpu used by bwa mem and sambamba (default: 8).'
    log.info '    --mem          INTEGER                 Size of memory used for alignment (in GB) (default: 32).'
    log.info '    --mem_sambamba          INTEGER                 Size of memory used by sambamba (in GB) (default: 1).'
    log.info '    --RG           STRING                  Samtools read group specification with "\t" between fields.'
    log.info '                                           e.g. --RG "PL:ILLUMINA\tDS:custom_read_group".'
    log.info '                                           Default: "PL:ILLUMINA".'
    log.info '    --fastq_ext      STRING                Extension of fastq files (default : fastq.gz)'
    log.info '    --bed        STRING                bed file with interval list'
    log.info '    --GATK_bundle        STRING                path to GATK bundle files (default : .)'
    log.info '    --GATK_folder        STRING                path to GATK GenomeAnalysisTK.jar file (default : .)'
    log.info '    --js        STRING                path to javascript interpreter k8'
    log.info '    --postaltjs        STRING                path to postalignment javascript bwa-postalt.js'
    log.info ""
    log.info "Flags:"
    log.info '--trim                    enable adapter sequence trimming'
    log.info '--indel_realignment                    performs local indel realignment (default: no).'
    log.info '--recalibration                    performs base quality score recalibration (GATK)'
    log.info '--alt                    enable alternative contig handling (for reference genome hg38)'
    log.info ''
    exit 0
}else {
  /* Software information */
  log.info "input_folder=${params.input_folder}"
  log.info "ref=${params.ref}"
  log.info "cpu=${params.cpu}"
  log.info "mem=${params.mem}"
  log.info "RG=${params.RG}"
  log.info "fastq_ext=${params.fastq_ext}"
  log.info "output_folder=${params.output_folder}"
  log.info "bed=${params.bed}"
  log.info "GATK_bundle=${params.GATK_bundle}"
  log.info "GATK_folder=${params.GATK_folder}"
  log.info "postaltjs=${params.postaltjs}"
  log.info "indel_realignment=${params.indel_realignment}"
  log.info "recalibration=${params.recalibration}"
  log.info "alt=${params.alt}"
  log.info "trim=${params.trim}"
  log.info "help=${params.help}"
}

//read files
ref = file(params.ref)
ref_fai = file( params.ref+'.fai' )
ref_sa = file( params.ref+'.sa' )
ref_bwt = file( params.ref+'.bwt' )
ref_ann = file( params.ref+'.ann' )
ref_amb = file( params.ref+'.amb' )
ref_pac = file( params.ref+'.pac' )
ref_alt = file( params.ref+'.alt' )
postaltjs = file( params.postaltjs )

import groovy.io.FileType
mode = 'fastq'
//inputsubfolders = file(params.input_folder).listFiles()
if(params.input_file){
	Channel.fromPath("${params.input_file}")
			.splitCsv()
			.map { row -> tuple("${row[0]}" , "${row[1]}" , file("${row[2]}"), file("${row[3]}")) }
			.set { readPairstmp }
	
	readPairs2merge = readPairstmp.groupTuple(by: 0)
	//bam_bai_files2merge.subscribe { row -> println "${row}" }
	single   = Channel.create()
	multiple = Channel.create()
	multiple1 = Channel.create()
	multiple2 = Channel.create()
	readPairs2merge.choice( single,multiple ) { a -> a[1].size() == 1 ? 0 : 1 }
	single2 = single.map { row -> tuple(row[0] , 1 , row[1][0], row[2][0], row[3][0])  }
      	//	.subscribe { row -> println "${row}" }
	multiple.separate(multiple1,multiple2){ row -> [ [row[0] , 2 ,  row[1][0], row[2][0], row[3][0]] , [row[0] , 2 , row[1][1], row[2][1],  row[3][1]] ] }
	readPairs=single2.concat(multiple1 ,multiple2 )
	//readPairs.subscribe { row -> println "${row}" }

	//Tests
	//Channel.fromPath("${params.input_file}")
        //                .splitCsv()
        //                .map { row -> tuple("${row[0]}" , "${row[1]}" , file("${row[2]}"), file("${row[3]}")) }
        //                .set {bam_bai_filesTest  }
	//bam_bai_files2merge = bam_bai_filesTest.groupTuple(by: 0)

	//queue1 = Channel.create()
	//queue2 = Channel.create()
	//bam_bai_files2merge.choice( queue1, queue2 ) { a -> a[1].size() == 1 ? 0 : 1 }
	//queue1.subscribe { row -> println "${row}" }
	//queue2.subscribe { row -> println "${row}" }

}else{
   if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.fastq_ext}/ }.size() > 0){
    	println "fastq files found, proceed with alignment"
	readPairs = Channel.fromFilePairs("${params.input_folder}/*_{1,2}*.${params.fastq_ext}")
			   .map { row -> tuple(row[0] , "" , row[1][0], row[1][1]) }
	Channel.fromFilePairs("${params.input_folder}/*{1,2}*.${params.fastq_ext}")
                           .map { row -> tuple(row[0] , "" , row[1][0], row[1][1]) }
			   .subscribe { row -> println "${row}" }
   }else{
    	if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0){
        	println "BAM files found, proceed with realignment"; mode ='bam'; files = Channel.fromPath( params.input_folder+'/*.bam' )
        }else{
        	println "ERROR: input folder contains no fastq nor BAM files"; System.exit(0)
        }
   }
}

//          readPairs = readPairs.concat( Channel.fromFilePairs("${folder}/*${params.fastq_ext}") ).subscribe { println it }

if(mode=='bam'){
    process bam_realignment {
        cpus params.cpu
        memory params.mem+'G'
        tag { file_tag }
        
        input:
        file infile from files
	file ref
	file ref_fai
	file ref_sa
	file ref_bwt
	file ref_ann
	file ref_amb
	file ref_pac
	file ref_alt
     
        output:
	set val(file_tag_new), file("${file_tag_new}.bam*") into bam_bai_files
	if( (params.recalibration==null)&(params.indel_realignment==null) ) publishDir params.output_folder, mode: 'move'

        shell:
	file_tag = infile.baseName
	file_tag_new=file_tag+'_realigned'
	if(params.trim) file_tag_new=file_tag_new+'_trimmed'
	if(params.alt)  file_tag_new=file_tag_new+'_alt'
	
	if(params.alt==null){
	  ignorealt='-j'
	  postalt=''
	}else{
	  ignorealt=''
	  postalt='k8 bwa-postalt.js '+params.ref+'.alt |'
	}
	if(params.trim==null){
	  preproc=''
	}else{
	  	
	  preproc='AdapterRemoval --interleaved --file1 /dev/stdin --output1 /dev/stdout |'
	}
	bwa_threads  = params.cpu.intdiv(2) - 1
	sort_threads = params.cpu.intdiv(2) - 1
	sort_mem     = params.mem.intdiv(4)
	'''
        set -o pipefail
        samtools collate -uOn 128 !{file_tag}.bam tmp_!{file_tag} | samtools fastq - | !{preproc} bwa mem !{ignorealt} -M -t!{bwa_threads} -R "@RG\\tID:!{file_tag}\\tSM:!{file_tag}\\t!{params.RG}" -p !{params.ref} - | !{postalt} samblaster --addMateTags | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t !{sort_threads} -m !{sort_mem}G --tmpdir=!{file_tag}_tmp -o !{file_tag_new}.bam /dev/stdin
        '''
    }
}
if(mode=='fastq'){
    println "fastq mode"
        
    process fastq_alignment {

        cpus params.cpu
        memory params.mem+'GB'    
        tag { "${file_tag}_${read_group}" }
        
        input:
        set val(file_tag), val(nb_rgs), val(read_group), file(pair1), file(pair2) from readPairs
	file ref
	file ref_fai
	file ref_sa
	file ref_bwt
	file ref_ann
	file ref_amb
	file ref_pac
	file ref_alt
	file postaltjs
                 
        output:
	set val(file_tag_new), val(nb_rgs), val(read_group),  file("${file_tag_new}*.bam*") into bam_bai_files
	if( (params.recalibration==null)&(params.indel_realignment==null) ) publishDir params.output_folder, mode: 'copy'

        shell:
	pair = [pair1,pair2]
	file_tag_new=file_tag 
	//+"_${read_group}"
	println file_tag_new
	bwa_threads  = params.cpu.intdiv(2) - 1
        sort_threads = params.cpu.intdiv(2) - 1
        sort_mem     = params.mem.intdiv(4)
	if(params.trim) file_tag_new=file_tag_new+'_trimmed'
        if(params.alt)  file_tag_new=file_tag_new+'_alt'	
	if(nb_rgs==1){
		file_tag_new=file_tag_new+"_${read_group}"
		compsort=" sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t ${sort_threads} -m ${sort_mem}G --tmpdir=${file_tag}_tmp -o ${file_tag_new}.bam /dev/stdin"
	}else{
		file_tag_new=file_tag_new
		compsort=" sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -n -t ${sort_threads} -m ${sort_mem}G --tmpdir=${file_tag}_tmp -o ${file_tag_new}_${read_group}.bam /dev/stdin"
	}
        if(params.alt==null){
          ignorealt='-j'
          postalt=''
        }else{
          ignorealt=''
          postalt='k8 bwa-postalt.js '+params.ref+'.alt |'
        }
	if(params.trim==null){
		'''
        	set -o pipefail
		bwa mem !{ignorealt} -M -t!{bwa_threads} -R "@RG\\tID:!{file_tag}.!{read_group}\\tSM:!{file_tag}\\t!{params.RG}" !{params.ref} !{pair[0]} !{pair[1]} | !{postalt} samblaster -M --addMateTags | !{compsort}
		'''
 	}else{
		'''
        	set -o pipefail
		AdapterRemoval --file1 !{pair[0]} --file2 !{pair[1]} --interleaved-output --output1 /dev/stdout | bwa mem !{ignorealt} -M -t!{bwa_threads} -R "@RG\\tID:!{file_tag}.!{read_group}\\tSM:!{file_tag}\\t!{params.RG}" -p !{params.ref} - | !{postalt} samblaster -M --addMateTags | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t !{sort_threads} -m !{sort_mem}G --tmpdir=!{file_tag}_tmp -o !{file_tag_new}.bam /dev/stdin
		'''
	}
	
     }
}

single_bam   = Channel.create()
multiple_bam0 = Channel.create()
bam_bai_files.choice( single_bam,multiple_bam0 ) { a -> a[1] == 1 ? 0 : 1 }

( mult2count, multiple_bam ) = multiple_bam0.into( 2 )

nmult = mult2count.count().println()
if(nmult>0 ){
	bam2merge = multiple_bam.groupTuple(by: 0)
			 .map { row -> tuple(row[0] , row[1][0] , row[2], row[3][0] , row[3][1] , null ,  null  ) }
}else{
	bam2merge = Channel.create()	
}
//( bam2merge, bam2mergeB ) = bam2merge1.into( 2 )

//bam2mergeB.subscribe { row -> println "${row}" }

process merge {
            cpus params.cpu
            memory params.mem+'G'
            tag { file_tag }
            input:
            set val(file_tag), val(nb_rgs), val(read_group),  file(bam1), file(bam2), file(bai1), file(bai2) from bam2merge
            output:
            set val(file_tag_new), file("${file_tag_new}.bam"), file("${file_tag_new}.bam.bai") into bam_bai_merged
             if( (params.recalibration==null)&(params.indel_realignment==null) ) publishDir params.output_folder, mode: 'move'

            shell:
            file_tag_new=file_tag+"_${read_group[0]}-${read_group[1]}_merged"
	    merge_threads  = params.cpu.intdiv(2) - 1
	    sort_threads = params.cpu.intdiv(2) - 1
            sort_mem     = params.mem.intdiv(4)
            '''
	    sambamba merge -t !{merge_threads} -l 0 /dev/stdout !{bam1} !{bam2} |  sambamba view -h /dev/stdin | samblaster -M --addMateTags | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t !{sort_threads} -m !{sort_mem}G --tmpdir=!{file_tag}_tmp -o !{file_tag_new}.bam /dev/stdin
            '''
}

bam_bai_files=single_bam.map { row -> tuple(row[0],row[3][0],row[3][1] ) }
			.concat(bam_bai_merged)
//                        .subscribe { row -> println "${row}" }


if(params.indel_realignment){ //Note: deprecated in gatk4
        // Local realignment around indels
        process indel_realignment {
            cpus params.cpu
            memory params.mem+'G'
            tag { file_tag }
            input:
	    set val(file_tag), file("${file_tag}.bam"), file("${file_tag}.bam.bai") from bam_bai_files
            output:
            file("${file_tag_new}_target_intervals.list") into indel_realign_target_files
            set val(file_tag_new), file("${file_tag_new}.bam"), file("${file_tag_new}.bam.bai") into bam_bai_files2
	    if(params.recalibration==null) publishDir params.output_folder, mode: 'move'
	    
            shell:
	    file_tag_new=file_tag+'_indelrealigned'
            '''
	    indelsvcf=`ls !{params.GATK_bundle}/*indels*.vcf* | grep -v ".tbi" | grep -v ".idx"`
	    knowncom=''
	    for ll in $indelsvcf; do knowncom=$knowncom' -known '$ll; done
            gatk RealignerTargetCreator --java-options "-Xmx!{params.mem}G" -R !{params.ref} -I !{file_tag}.bam $knowncom -o !{file_tag_new}_target_intervals.list
            gatk IndelRealigner --java-options "-Xmx!{params.mem}G" -R !{params.ref} -I !{file_tag}.bam -targetIntervals !{file_tag_new}_target_intervals.list $knowncom -o !{file_tag_new}.bam
	    mv !{file_tag_new}.bai !{file_tag_new}.bam.bai
            '''
        }
}else{
    if(params.recalibration){
        bam_bai_files2 = bam_bai_files
    }
}

if(params.recalibration){
// base quality score recalibration
   process base_quality_score_recalibration {
    cpus params.cpu
    memory params.mem+'G'
    tag { file_tag }
        
    input:
    set val(file_tag), file("${file_tag}.bam"), file("${file_tag}.bam.bai") from bam_bai_files2
    output:
    file("*_recal.table") into recal_table_files
    file("*plots.pdf") into recal_plots_files
    set val(file_tag_new), file("${file_tag_new}.bam"), file("${file_tag_new}.bam.bai") into recal_bam_bai_files
    publishDir params.output_folder, mode: 'move'

    shell:
    file_tag_new=file_tag+'_BQSRecalibrated'
    '''
    indelsvcf=(`ls !{params.GATK_bundle}/*indels*.vcf* | grep -v ".tbi" | grep -v ".idx"`)
    dbsnpvcfs=(`ls !{params.GATK_bundle}/*dbsnp*.vcf* | grep -v ".tbi" | grep -v ".idx"`)
    dbsnpvcf=${dbsnpvcfs[@]:(-1)}
    knownSitescom=''
    for ll in $indelsvcf; do knownSitescom=$knownSitescom' --known-sites '$ll; done
    knownSitescom=$knownSitescom' --known-sites '$dbsnpvcf
    gatk BaseRecalibrator --java-options "-Xmx!{params.mem}G" -R !{params.ref} -I !{file_tag}.bam $knownSitescom -O !{file_tag}_recal.table
    gatk ApplyBQSR --java-options "-Xmx!{params.mem}G" -R !{params.ref} -I !{file_tag}.bam --bqsr-recal-file !{file_tag}_recal.table -O !{file_tag_new}.bam
    gatk BaseRecalibrator --java-options "-Xmx!{params.mem}G" -R !{params.ref} -I !{file_tag}.bam $knownSitescom -O !{file_tag_new}_recal.table		
    gatk AnalyzeCovariates --java-options "-Xmx!{params.mem}G" -before !{file_tag}_recal.table -after !{file_tag_new}_recal.table -plots !{file_tag_new}_recalibration_plots.pdf	
    mv !{file_tag_new}.bai !{file_tag_new}.bam.bai
    '''
    }
}
