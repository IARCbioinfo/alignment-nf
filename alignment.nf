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
params.ref          = 'hg19.fasta'
params.cpu          = 8
params.mem          = 32
params.RG           = "PL:ILLUMINA"
params.fastq_ext    = "fastq.gz"
params.suffix1      = "_1"
params.suffix2      = "_2"
params.output_folder = "."
params.bed          = ""
params.GATK_bundle  = "bundle"
params.GATK_folder  = "."
params.js           = "k8"
params.postaltjs    = "bwa-postalt.js"
params.indel_realignment = null
params.recalibration = null
params.help         = null
params.alt          = null
params.trim         = null


log.info ""
log.info "--------------------------------------------------------"
log.info "  alignment-nf 1.0.0: alignment/realignment workflow for whole exome/whole genome sequencing "
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
    log.info '    --suffix1        STRING                Suffix of fastq files 1 (default : _1)'
    log.info '    --suffix2        STRING                Suffix of fastq files 2 (default : _2)'
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
  log.info "suffix1=${params.suffix1}"
  log.info "suffix2=${params.suffix2}"
  log.info "output_folder=${params.output_folder}"
  log.info "bed=${params.bed}"
  log.info "GATK_bundle=${params.GATK_bundle}"
  log.info "GATK_folder=${params.GATK_folder}"
  log.info "js=${params.js}"
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
	file ref
	file ref_fai
	file ref_sa
	file ref_bwt
	file ref_ann
	file ref_amb
	file ref_pac
	file ref_alt
     
        output:
	set val(file_tag_new), file("${file_tag_new}.bam") into bam_files
	file("${file_tag_new}.bam.bai") into bai_files
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
	  postalt=params.js+' '+params.postaltjs+' '+params.ref+'.alt |'
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
	file ref
	file ref_fai
	file ref_sa
	file ref_bwt
	file ref_ann
	file ref_amb
	file ref_pac
	file ref_alt
                 
        output:
	set val(file_tag_new), file("${file_tag_new}.bam") into bam_files
	file("${file_tag_new}.bam.bai") into bai_files
	if( (params.recalibration==null)&(params.indel_realignment==null) ) publishDir params.output_folder, mode: 'move'

        shell:
        file_tag = pair[0].name.replace("${params.suffix1}.${params.fastq_ext}","")
	file_tag_new=file_tag
	if(params.trim) file_tag_new=file_tag_new+'_trimmed'
	if(params.alt)  file_tag_new=file_tag_new+'_alt'
	if(params.alt==null){
	  ignorealt='-j'
	  postalt=''
	}else{
	  ignorealt=''
	  postalt=params.js+' '+params.postaltjs+' '+params.ref+'.alt |'
	}
	bwa_threads  = params.cpu.intdiv(2) - 1
	sort_threads = params.cpu.intdiv(2) - 1
	sort_mem     = params.mem.intdiv(4)
	if(params.trim==null){
		'''
        	set -o pipefail
		bwa mem !{ignorealt} -M -t!{bwa_threads} -R "@RG\\tID:!{file_tag}\\tSM:!{file_tag}\\t!{params.RG}" !{params.ref} !{pair[0]} !{pair[1]} | !{postalt} samblaster --addMateTags | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t !{sort_threads} -m !{sort_mem}G --tmpdir=!{file_tag}_tmp -o !{file_tag_new}.bam /dev/stdin
		'''
 	}else{
		'''
        	set -o pipefail
		AdapterRemoval --file1 !{pair[0]} --file2 !{pair[1]} --interleaved-output --output1 /dev/stdout | bwa mem !{ignorealt} -M -t!{bwa_threads} -R "@RG\\tID:!{file_tag}\\tSM:!{file_tag}\\t!{params.RG}" -p !{params.ref} - | !{postalt} samblaster --addMateTags | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t !{sort_threads} -m !{sort_mem}G --tmpdir=!{file_tag}_tmp -o !{file_tag_new}.bam /dev/stdin
		'''
	}
	
     }
}

if(params.indel_realignment){
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
	    if(params.recalibration==null) publishDir params.output_folder, mode: 'move'
	    
            shell:
	    file_tag_new=file_tag+'_indelrealigned'
            '''
	    indelsvcf=`ls !{params.GATK_bundle}/*indels*.vcf* | grep -v ".tbi" | grep -v ".idx"`
	    knowncom=''
	    for ll in $indelsvcf; do knowncom=$knowncom' -known '$ll; done
            java -Xmx!{params.mem}g -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt !{params.cpu} -R !{params.ref} -I !{file_tag}.bam $knowncom -o !{file_tag_new}_target_intervals.list
            java -Xmx!{params.mem}g -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T IndelRealigner -R !{params.ref} -I !{file_tag}.bam -targetIntervals !{file_tag_new}_target_intervals.list $knowncom -o !{file_tag_new}.bam
	    mv !{file_tag_new}.bai !{file_tag_new}.bam.bai
            '''
        }
}else{
    if(params.recalibration){
        bam_files2 = bam_files
	bai_files2 = bai_files
    }
}

if(params.recalibration){
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
    publishDir params.output_folder, mode: 'move'

    shell:
    file_tag_new=file_tag+'_BQSrecalibrated'
    '''
    indelsvcf=(`ls !{params.GATK_bundle}/*indels*.vcf* | grep -v ".tbi" | grep -v ".idx"`)
    dbsnpvcfs=(`ls !{params.GATK_bundle}/*dbsnp*.vcf* | grep -v ".tbi" | grep -v ".idx"`)
    dbsnpvcf=${dbsnpvcfs[@]:(-1)}
    knownSitescom=''
    for ll in $indelsvcf; do knownSitescom=$knownSitescom' -knownSites '$ll; done
    knownSitescom=$knownSitescom' -knownSites '$dbsnpvcf
    java -Xmx!{params.mem}g -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T BaseRecalibrator -nct !{params.cpu} -R !{params.ref} -I !{file_tag}.bam $knownSitescom -L !{params.bed} -o !{file_tag_new}_recal.table
    java -Xmx!{params.mem}g -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T BaseRecalibrator -nct !{params.cpu} -R !{params.ref} -I !{file_tag}.bam $knownSitescom -BQSR !{file_tag_new}_recal.table -L !{params.bed} -o !{file_tag_new}_post_recal.table		
    java -Xmx!{params.mem}g -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T AnalyzeCovariates -R !{params.ref} -before !{file_tag_new}_recal.table -after !{file_tag_new}_post_recal.table -plots !{file_tag_new}_recalibration_plots.pdf	
    java -Xmx!{params.mem}g -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T PrintReads -nct !{params.cpu} -R !{params.ref} -I !{file_tag}.bam -BQSR !{file_tag_new}_recal.table -L !{params.bed} -o !{file_tag_new}.bam
    mv !{file_tag_new}.bai !{file_tag_new}.bam.bai
    '''
    }
}
