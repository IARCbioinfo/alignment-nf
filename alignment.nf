#! /usr/bin/env nextflow

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

params.input_folder = null
params.input_file   = null
params.ref          = 'hg19.fasta'
params.cpu          = 8
params.mem          = 32
params.RG           = "PL:ILLUMINA"
params.fastq_ext    = "fastq.gz"
params.suffix1      = "_1"
params.suffix2      = "_2"
params.output_folder = "."
params.bed          = ""
params.snp_vcf      = "dbsnp.vcf"
params.indel_vcf    = "Mills_1000G_indels.vcf"
params.postaltjs    = "bwa-postalt.js"
params.feature_file = 'NO_FILE'
params.mem_BQSR     = 10
params.cpu_BQSR     = 2
params.bwa_option_M  = null
params.recalibration = null
params.help         = null
params.alt          = null
params.trim         = null


log.info ""
log.info "--------------------------------------------------------"
log.info "  alignment-nf 1.1.0: alignment/realignment workflow for whole exome/whole genome sequencing "
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
    log.info '--input_folder   FOLDER              Folder containing BAM or fastq files to be aligned.'
    log.info '--ref            FILE                Reference fasta file (with index).'
    log.info ""
    log.info 'Optional arguments:'
    log.info '--input_file     STRING              Input file (comma-separated) with 3 columns:'
    log.info '                                     sample name, read_group_ID, and file path.'
    log.info '--output_folder  STRING              Output folder (default: .).'
    log.info '--cpu            INTEGER             Number of cpu used by bwa mem and sambamba (default: 8).'
    log.info '--mem            INTEGER             Size of memory used for alignment (in GB) (default: 32).'
    log.info '--RG             STRING              Samtools read group specification with "\t" between fields.'
    log.info '                                         e.g. --RG "PL:ILLUMINA\tDS:custom_read_group".'
    log.info '                                         Default: "PL:ILLUMINA".'
    log.info '--fastq_ext      STRING              Extension of fastq files (default: fastq.gz)'
    log.info '--suffix1        STRING              Suffix of fastq files 1 (default : _1)'
    log.info '--suffix2        STRING              Suffix of fastq files 2 (default : _2)'
    log.info '--bed            STRING              Bed file with interval list'
    log.info '--snp_vcf        STRING              Path to SNP VCF from GATK bundle (default: dbsnp.vcf)'
    log.info '--indel_vcf      STRING              Path to indel VCF from GATK bundle (default: Mills_1000G_indels.vcf)'
    log.info '--postaltjs      STRING              Path to postalignment javascript bwa-postalt.js'
    log.info '--feature_file   STRING              Path to feature file for qualimap (default: NO_FILE)'
    log.info '--mem_BQSR       INTEGER             Size of memory used for GATK BQSR (in GB) (default: 10)'
    log.info '--cpu_BQSR       INTEGER             Number of cpu used by GATK BQSR (default: 2)'
    log.info ""
    log.info "Flags:"
    log.info '--trim                               Enable adapter sequence trimming'
    log.info '--recalibration                      Performs base quality score recalibration (GATK)'
    log.info '--alt                                Enable alternative contig handling (for reference genome hg38)'
    log.info '--bwa_option_M                       Trigger the -M option in bwa and the corresponding compatibility option in samblaster'
    log.info ''
    exit 0
}else {
  /* Software information */
  log.info "input_folder=${params.input_folder}"
  log.info "input_file=${params.input_file}"
  log.info "ref=${params.ref}"
  log.info "cpu=${params.cpu}"
  log.info "mem=${params.mem}"
  log.info "RG=${params.RG}"
  log.info "fastq_ext=${params.fastq_ext}"
  log.info "suffix1= ${params.suffix1}"
  log.info "suffix2= ${params.suffix2}"
  log.info "output_folder=${params.output_folder}"
  log.info "bed=${params.bed}"
  log.info "snp_vcf=${params.snp_vcf}"
  log.info "indel_vcf=${params.indel_vcf}"
  log.info "postaltjs=${params.postaltjs}"
  log.info "feature_file=${params.feature_file}"
  log.info "mem_BQSR=${params.mem_BQSR}"
  log.info "cpu_BQSR=${params.cpu_BQSR}"
  log.info "recalibration=${params.recalibration}"
  log.info "alt=${params.alt}"
  log.info "trim=${params.trim}"
  log.info "bwa_option_M=${params.bwa_option_M}"
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
ref_dict= file( params.ref.replaceFirst(/fasta/, "").replaceFirst(/fa/, "") +'dict')
postaltjs = file( params.postaltjs )

//get know site VCFs from GATK bundle
known_snps         = file( params.snp_vcf )
known_snps_index   = file( params.snp_vcf+'.tbi' )
known_indels       = file( params.indel_vcf )
known_indels_index = file( params.indel_vcf+'.tbi' )

//qualimap feature file
qualimap_ff = file(params.feature_file)

mode = 'fastq'
if(params.input_file){
	Channel.fromPath("${params.input_file}")
			.splitCsv()
			.map { row -> tuple("${row[0]}" , "${row[1]}" , file("${row[2]}"), file("${row[3]}")) }
			.set { readPairstmp }
	
	readPairs2merge = readPairstmp.groupTuple(by: 0)
	single   = Channel.create()
	multiple = Channel.create()
	multiple1 = Channel.create()
	multiple2 = Channel.create()
	readPairs2merge.choice( single,multiple ) { a -> a[1].size() == 1 ? 0 : 1 }
	single2 = single.map { row -> tuple(row[0] , 1 , row[1][0], row[2][0], row[3][0])  }
	multiple.separate(multiple1,multiple2){ row -> [ [row[0] , 2 ,  row[1][0], row[2][0], row[3][0]] , [row[0] , 2 , row[1][1], row[2][1],  row[3][1]] ] }
	readPairs=single2.concat(multiple1 ,multiple2 )
}else{
   if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.fastq_ext}/ }.size() > 0){
    	println "fastq files found, proceed with alignment"
	readPairs = Channel.fromFilePairs(params.input_folder +"/*{${params.suffix1},${params.suffix2}}" +'.'+ params.fastq_ext)
			   .map { row -> [ row[0] , 1 , row[0] , row[1][0], row[1][1] ] }
   }else{
    	if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0){
        	println "BAM files found, proceed with realignment"; mode ='bam'; files = Channel.fromPath( params.input_folder+'/*.bam' )
        }else{
        	println "ERROR: input folder contains no fastq nor BAM files"; System.exit(0)
        }
   }
}

if(mode=='bam'){
    process bam_realignment {
        cpus params.cpu
        memory params.mem+'G'
        tag { file_tag }
        
        if(!params.recalibration) publishDir "${params.output_folder}/BAM/", mode: 'copy'

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
	file postaltjs
     
        output:
	set val(file_tag_new), val(1), val("RG"), file("${file_tag_new}*.bam*")  into bam_bai_files0

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
	if(params.bwa_option_M==null){
	  bwa_opt=''
    	  samblaster_opt=''
	}else{
	   bwa_opt='-M '
 	   samblaster_opt='-M '
	}
	bwa_threads  = [params.cpu.intdiv(2) - 1,1].max()
	sort_threads = [params.cpu.intdiv(2) - 1,1].max()
	sort_mem     = params.mem.intdiv(4)
	'''
        set -o pipefail
        samtools collate -uOn 128 !{file_tag}.bam tmp_!{file_tag} | samtools fastq - | !{preproc} bwa mem !{ignorealt} !{bwa_opt} -t!{bwa_threads} -R "@RG\\tID:!{file_tag}\\tSM:!{file_tag}\\t!{params.RG}" -p !{ref} - | !{postalt} samblaster !{samblaster_opt} --addMateTags --ignoreUnmated | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t !{sort_threads} -m !{sort_mem}G --tmpdir=!{file_tag}_tmp -o !{file_tag_new}.bam /dev/stdin
        '''
    }
}
if(mode=='fastq'){
    println "fastq mode"
        
    process fastq_alignment {
        cpus params.cpu
        memory params.mem+'GB'    
        tag { "${file_tag}_${read_group}" }
 	
	if(!params.recalibration){ publishDir "${params.output_folder}/BAM/", mode: 'copy', 
            saveAs: {filename -> 
                if (nb_rgs == 1) "$filename"
                else null
		}
	}

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
	set val(file_tag_new), val(nb_rgs), val(read_group),  file("${file_tag_new}*.bam*") into bam_bai_files0

        shell:
	pair = [pair1,pair2]
	file_tag_new=file_tag 
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
	if(params.bwa_option_M==null){
          bwa_opt=''
          samblaster_opt=''
        }else{
           bwa_opt='-M '
           samblaster_opt='-M '
        }
	if(params.trim==null){
		'''
        	set -o pipefail
		bwa mem !{ignorealt} -M -t!{bwa_threads} -R "@RG\\tID:!{file_tag}.!{read_group}\\tSM:!{file_tag}\\t!{params.RG}" !{ref} !{pair[0]} !{pair[1]} | !{postalt} samblaster -M --addMateTags | !{compsort}
		'''
 	}else{
		'''
        	set -o pipefail
		AdapterRemoval --file1 !{pair[0]} --file2 !{pair[1]} --interleaved-output --output1 /dev/stdout | bwa mem !{ignorealt} !{bwa_opt} -t!{bwa_threads} -R "@RG\\tID:!{file_tag}.!{read_group}\\tSM:!{file_tag}\\t!{params.RG}" -p !{ref} - | !{postalt} samblaster !{samblaster_opt} --addMateTags | !{compsort}
		'''
	}
     }
}

if(params.input_file){
	println "Merge"
	single_bam   = Channel.create()
	multiple_bam0 = Channel.create()
	bam_bai_files0.choice( single_bam,multiple_bam0 ) { a -> a[1] == 1 ? 0 : 1 }
	( mult2count, mult2QC, multiple_bam ) = multiple_bam0.into( 3 )
	
	//QC on each run
	process qualimap_multi {
	    cpus params.cpu
	    memory params.mem+'G'
	    tag { "${file_tag}_${read_group}" }

	    publishDir "${params.output_folder}/QC/BAM/qualimap/", mode: 'copy'

	    input:
	    set val(file_tag), val(nb_rgs), val(read_group),  file(bambai) from mult2QC
	    file qff from qualimap_ff

	    output:
	    file ("${file_tag}_${read_group}") into qualimap_multi_results
	    file ("${file_tag}_${read_group}.stats.txt") into flagstat_multi_results

	    shell:
	    bam=bambai[0]
	    feature = qff.name != 'NO_FILE' ? "--feature-file $qff" : ''
	    '''
	    sambamba sort -t !{params.cpu} -m !{params.mem}G --tmpdir=!{file_tag}_tmp -o !{file_tag}_!{read_group}_COsorted.bam !{bam}
	    qualimap bamqc -nt !{params.cpu} !{feature} --skip-duplicated -bam !{file_tag}_!{read_group}_COsorted.bam --java-mem-size=!{params.mem}G -outdir !{file_tag}_!{read_group} -outformat html
	    sambamba flagstat -t !{params.cpu} !{bam} > !{file_tag}_!{read_group}.stats.txt
	    '''
	}

	process multiqc_multi {
	    cpus 2
	    memory '1G'

	    publishDir "${params.output_folder}/QC/BAM/qualimap", mode: 'copy'

	    input:
	    file qualimap_results from qualimap_multi_results.collect()
	    file flagstat_results from flagstat_multi_results.collect()

	    output:
	    file("*report.html") into multi_output
	    file("multiqc_multiplex_qualimap_flagstat_report_data/") into multi_output_data

	    shell:
	    '''
	    multiqc . -n multiqc_multiplex_qualimap_flagstat_report.html
	    '''
	}

	//merge runs
	//nmult = mult2count.toList().size() //count().subscribe{ print "$it" }
	//println nmult
	//if( nmult >0 ){
		//println "BAMs from multiple runs detected"
		bam2merge = multiple_bam.groupTuple(by: 0)
			 .map { row -> tuple(row[0] , row[1][0] , row[2], row[3][0] , row[3][1] , null ,  null  ) }
	//}else{
	//	println "No BAMs from multiple runs detected"
	//	bam2merge = Channel.create()	
	//}

	process merge {
            cpus params.cpu
            memory params.mem+'G'
            tag { file_tag }
            if(!params.recalibration) publishDir "$params.output_folder/BAM/", mode: 'copy', pattern: "*.bam*"

            input:
            set val(file_tag), val(nb_rgs), val(read_group),  file(bam1), file(bam2), file(bai1), file(bai2) from bam2merge

            output:
            set val(file_tag_new), file("${file_tag_new}.bam"), file("${file_tag_new}.bam.bai") into bam_bai_merged

            shell:
            file_tag_new=file_tag+"_${read_group[0]}-${read_group[1]}_merged"
	    merge_threads  = params.cpu.intdiv(2) - 1
	    sort_threads = params.cpu.intdiv(2) - 1
            sort_mem     = params.mem.intdiv(2)
            '''
	    sambamba merge -t !{merge_threads} -l 0 /dev/stdout !{bam1} !{bam2} |  sambamba view -h /dev/stdin | samblaster -M --addMateTags | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t !{sort_threads} -m !{sort_mem}G --tmpdir=!{file_tag}_tmp -o !{file_tag_new}.bam /dev/stdin
            '''
	}

	bam_bai_files=single_bam.map { row -> tuple(row[0],row[3][0],row[3][1] ) }
			.concat(bam_bai_merged)
}else{
	bam_bai_files = bam_bai_files0.map { row -> tuple(row[0],row[3][0],row[3][1] ) }
}

if(params.recalibration){
println "BQSR"

// base quality score recalibration
   process base_quality_score_recalibration {
    cpus params.cpu_BQSR
    memory params.mem_BQSR+'G'
    tag { file_tag }
    
    publishDir "$params.output_folder/BAM/", mode: 'copy', pattern: "*bam*"
    publishDir "$params.output_folder/QC/BAM/BQSR/", mode: 'copy',
	saveAs: {filename -> 
		if (filename.indexOf("table") > 0) "$filename"
		else if (filename.indexOf("plots") > 0) "$filename"
		else null
	}

    input:
    set val(file_tag), file("${file_tag}.bam"), file("${file_tag}.bam.bai") from bam_bai_files
    file known_snps
    file known_snps_index
    file known_indels
    file known_indels_index
    file ref
    file ref_fai
    file ref_dict

    output:
    file("*_recal.table") into recal_table_files
    file("*plots.pdf") into recal_plots_files
    set val(file_tag_new), file("${file_tag_new}.bam"), file("${file_tag_new}.bam.bai") into final_bam_bai_files

    shell:
    file_tag_new=file_tag+'_BQSRecalibrated'
    '''
    gatk BaseRecalibrator --java-options "-Xmx!{params.mem_BQSR}G" -R !{ref} -I !{file_tag}.bam --known-sites !{known_snps} --known-sites !{known_indels} -O !{file_tag}_recal.table
    gatk ApplyBQSR --java-options "-Xmx!{params.mem_BQSR}G" -R !{ref} -I !{file_tag}.bam --bqsr-recal-file !{file_tag}_recal.table -O !{file_tag_new}.bam
    gatk BaseRecalibrator --java-options "-Xmx!{params.mem_BQSR}G" -R !{ref} -I !{file_tag_new}.bam --known-sites !{known_snps} --known-sites !{known_indels} -O !{file_tag_new}_recal.table		
    gatk AnalyzeCovariates --java-options "-Xmx!{params.mem_BQSR}G" -before !{file_tag}_recal.table -after !{file_tag_new}_recal.table -plots !{file_tag_new}_recalibration_plots.pdf	
    mv !{file_tag_new}.bai !{file_tag_new}.bam.bai
    '''
    }
}else{
	final_bam_bai_files = bam_bai_files
	recal_table_files = Channel.from ( 'NOFILE1', 'NOFILE2' )
}


process qualimap_final {
    cpus params.cpu
    memory params.mem+'G'
    tag { file_tag }

    publishDir "${params.output_folder}/QC/BAM/qualimap/", mode: 'copy'

    input:
    set val(file_tag), file(bam), file(bai) from final_bam_bai_files
    file qff from qualimap_ff

    output:
    file ("${file_tag}") into qualimap_results
    file ("${file_tag}.stats.txt") into flagstat_results

    shell:
    feature = qff.name != 'NO_FILE' ? "--feature-file $qff" : ''
    '''
    qualimap bamqc -nt !{params.cpu} !{feature} --skip-duplicated -bam !{bam} --java-mem-size=!{params.mem}G -outdir !{file_tag} -outformat html
    sambamba flagstat -t !{params.cpu} !{bam} > !{file_tag}.stats.txt
    '''
}

process multiqc_final {
    cpus 2
    memory '2G'

    publishDir "${params.output_folder}/QC/BAM/", mode: 'copy'

    input:
    file qualimap_results from qualimap_results.collect()
    file flagstat_results from flagstat_results.collect()
    file BQSR_results from recal_table_files.collect()

    output:
    file("*report.html") into final_output
    file("multiqc_qualimap_flagstat_BQSR_report_data/") into final_output_data

    shell:
    '''
    multiqc . -n multiqc_qualimap_flagstat_BQSR_report.html
    '''
}

// Display completion message
workflow.onComplete {
  log.info "N E X T F L O W  ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
  //log.info "iarcbioinfo/alignment-nf ~ " + this.grabRevision() + (workflow.commitId ? " [${workflow.commitId}]" : "")
  log.info "Completed at: " + workflow.complete
  log.info "Duration    : " + workflow.duration
  log.info "Success     : " + workflow.success
  log.info "Exit status : " + workflow.exitStatus
  log.info "Error report: " + (workflow.errorReport ?: '-')
}
