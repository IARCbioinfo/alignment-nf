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
//params.output_folder = "./results"
params.bed          = ""
params.snp_vcf      = "dbsnp.vcf"
params.indel_vcf    = "Mills_1000G_indels.vcf"
params.postaltjs    = "NO_FILE"
params.feature_file = 'NO_FILE'
params.mem_BQSR     = 10
params.cpu_BQSR     = 2
params.multiqc_config = 'NO_FILE'
params.adapterremoval_opt = ""
params.bwa_mem      = "bwa-mem2 mem"
params.bwa_option_M  = null
params.recalibration = null
params.help         = null
params.alt          = null
params.trim         = null
//bwakit directory
params.bwakit_root      = '/opt/conda/envs/alignment-nf/share/bwakit-0.7.15-1/'

//new variables
params.output_type         = "cram" //default output type is cram
params.cram_ref = "def.cram.fasta"
params.bazam = "bazam" //save location in the docker or singularity container

//we display the header of the tool
log.info IARC_Header()
log.info tool_header()

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
    log.info '--input_file     STRING              Input file (comma-separated) with 4 columns:'
    log.info '                                     SM(sample name), RG (read_group_ID), pair1 (path to fastq pair 1), '
    log.info '                                     and pair2 (path to fastq pair 2).'
    log.info '--cram_ref       STRING              Path to CRAM reference in fasta format to perform realigment, the reference must be indexed (samtools faidx)'
    log.info '--output_folder  STRING              Output folder (default: "./results").'
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
    log.info '--bwakit_root                        Root directory of bwakit'
    log.info '--feature_file   STRING              Path to feature file for qualimap (default: NO_FILE)'
    log.info '--mem_BQSR       INTEGER             Size of memory used for GATK BQSR (in GB) (default: 10)'
    log.info '--cpu_BQSR       INTEGER             Number of cpu used by GATK BQSR (default: 2)'
    log.info '--multiqc_config STRING              Config yaml file for multiqc (default : none)'
    log.info '--adapterremoval_opt STRING          Command line options for AdapterRemoval (default : none)'
    log.info '--bwa_mem        STRING              bwa-mem command (default: "bwa-mem2 mem", alternative is "bwa mem")'
    log.info '--bazam          STRING              Path to bazam program (default: "bazam")'
    log.info ""
    log.info "Flags:"
    log.info '--trim                               Enable adapter sequence trimming'
    log.info '--recalibration                      Performs base quality score recalibration (GATK)'
    log.info '--alt                                Enable alternative contig handling (for reference genome hg38)'
    log.info '--bwa_option_M                       Trigger the -M option in bwa and the corresponding compatibility option in samblaster'
    log.info '--output_type  STRING                Default is CRAM, but BAM is still supported.'
    log.info ''
    exit 0
}else {
  /* Software information */
//we print the parameters
log.info "\n"
log.info "-\033[2m------------------Calling PARAMETERS--------------------\033[0m-"
log.info params.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
log.info "-\033[2m--------------------------------------------------------\033[0m-"
log.info "\n"
//todo: add software versions
}

//multiqc config file
ch_config_for_multiqc = file(params.multiqc_config)

//read files
ref     = file(params.ref)
ref_fai = file(params.ref+'.fai')
ref_sa  = file(params.ref+'.sa')
ref_bwt =  file(params.ref+'.bwt')
ref_ann =  file(params.ref+'.ann')
ref_amb =  file(params.ref+'.amb')
ref_pac =  file(params.ref+'.pac')
ref_dict=  file(params.ref.replaceFirst(/fasta/, "").replaceFirst(/fa/, "") +'dict')

ref_cram = file(params.cram_ref)

if(params.bwa_mem!="bwa-mem2 mem"){
  ref_0123 = file('NO_0123')
  ref_bwt8bit = file('NO_bwt8bit')
}else{
  ref_0123 = file(params.ref+'.0123')
  //ref_bwt2bit = params.ref+'.bwt.2bit.64'
  //ref_bwt8bit = file(params.ref+'.bwt.8bit.32')
  ref_bwt8bit = file(params.ref+'.bwt.2bit.64')
}
//bwa-mem2 files
if(params.alt){
  ref_alt = file(params.ref+'.alt')
}else{
  ref_alt = file('NO_ALT')
}
postaltjs = file( params.postaltjs )

//get know site VCFs from GATK bundle
known_snps         = file( params.snp_vcf )
known_snps_index   = file( params.snp_vcf+'.tbi' )
known_indels       = file( params.indel_vcf )
known_indels_index = file( params.indel_vcf+'.tbi' )

//qualimap feature file
qualimap_ff = file(params.feature_file)

assert (params.input_file != null | params.input_folder != null) : "please specify input_file or input_folder"
//we create a channel_for_reference
ch_ref=Channel.value(file(params.ref)).ifEmpty{exit 1, "reference file not found: ${params.ref}"}

mode = 'fastq'
if(params.input_file){
Channel.fromPath("${params.input_file}")
			.splitCsv(header: true, sep: '\t', strip: true)
			.map { row -> [row.SM , "_"+row.RG , file(row.pair1), file(row.pair2) ] }
      .into{readPairs0;readPairs4group}

  readPairsgrouped = readPairs4group.groupTuple(by: 0)
	                                 .map{ a -> [a[0],a[1].size(),a[1],a[2],a[3]] }

	readPairs = readPairsgrouped.map{ a -> [a[0],a[1]] }
                            .cross( readPairs0 )
                            .map{a -> [a[1][0],a[0][1],a[1][1],a[1][2],a[1][3] ] }
}else{
   if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.fastq_ext}/ }.size() > 0){
    	println "fastq files found, proceed with alignment"
	readPairs = Channel.fromFilePairs(params.input_folder +"/*{${params.suffix1},${params.suffix2}}" +'.'+ params.fastq_ext)
			   .map { row -> [ row[0] , 1 , "" , row[1][0], row[1][1] ] }
   }else{
      //we try BAM files
    	if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0){
        	println "BAM files found, proceed with realignment";
          mode ='bam'
          bams = Channel.fromPath( params.input_folder+'/*.bam' )
                  .map {path -> [ path.name.replace(".bam",""),path]}
          bams_index = Channel.fromPath( params.input_folder+'/*.bam.bai')
                  .map {  path -> [ path.name.replace(".bam.bai",""), path ] }
          //we create the chanel
          files = bams.join(bams_index)

        }else{
          //we try CRAM files
           if(file(params.input_folder).listFiles().findAll { it.name ==~ /.*cram/ }.size() > 0){
               println "CRAM files found, proceed with realignment";
               mode ='cram';
               crams = Channel.fromPath( params.input_folder+'/*.cram')
                        .map {path -> [ path.name.replace(".cram",""),path]}
               crams_index = Channel.fromPath( params.input_folder+'/*.cram.crai')
                            .map {  path -> [ path.name.replace(".cram.crai",""), path ] }
                //we create the chanel prefix .cram .cram.crai
                files = crams.join(crams_index)
                 if(params.cram_ref == "def.cram.fasta"){
                   println "We detected CRAM files for realignment, but the CRAM reference was not set (--cram_ref)"; System.exit(1)
               }

         }else{
        	println "ERROR: input folder contains no fastq nor BAM/CRAM files"; System.exit(1)
         }
      }
   }
}

if(mode=='bam' || mode=='cram'){
  process bam_realignment {
    tag { file_tag }
    label 'realn_cram_bam'
    //if(!params.recalibration) publishDir "${params.output_folder}/BAM/", mode: 'copy'

	input:
  //file infile, index from files
  set val(id), file(infile), file(infile_index) from files
	file ref
  file ref_sa
  file ref_bwt
  file ref_ann
  file ref_amb
  file ref_pac
  file ref_0123
  file ref_bwt8bit
  file ref_alt
	file postaltjs
  file ref_cram

  output:
	set val(file_tag), file("${file_tag_new}*.bam"), file("${file_tag_new}*.bai")  into bam_bai_files, bam_bai_to_cram_files

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
	  //postalt='k8 bwa-postalt.js '+ref+'.alt |'
    //heng li code $root/k8 $root/bwa-postalt.js $hla_pre$ARGV[0].alt
    postalt="${params.bwakit_root}/k8 ${params.bwakit_root}/bwa-postalt.js "+ref+".alt |"
	}
	if(params.trim==null){
	  preproc=''
	}else{
	  preproc="AdapterRemoval ${params.adapterremoval_opt} --interleaved --file1 /dev/stdin --output1 /dev/stdout |"
	}
	if(params.bwa_option_M==null){
	  bwa_opt=''
    samblaster_opt=''
	}else{
	   bwa_opt='-M '
 	   samblaster_opt='-M '
	}
	//bwa_threads  = [params.cpu.intdiv(2) - 1,1].max()
//	sort_threads = [params.cpu.intdiv(2) - 1,1].max()
	sort_mem     = params.mem.div(4)
  bwa_threads = task.cpus
  sort_threads = task.cpus


  if(mode == 'bam'){
    //todo: kept samtools collect option [def:bazam] for CRAM/BAM
  	'''
    set -o pipefail
    bazam -n 1 -Xms2G -Xmx20G  -bam !{file_tag}.bam | \\
    !{preproc} !{params.bwa_mem} !{ignorealt} !{bwa_opt} -t!{bwa_threads} -R "@RG\\tID:!{file_tag}\\tSM:!{file_tag}\\t!{params.RG}" -p !{ref} - | \\
    !{postalt} samblaster !{samblaster_opt} --addMateTags --ignoreUnmated | \\
    sambamba view -S -f bam -l 0 /dev/stdin | \\
    sambamba sort -t !{sort_threads} -m !{sort_mem}G --tmpdir=!{file_tag}_tmp -o !{file_tag_new}.bam /dev/stdin
    '''
  }else{
    '''
    set -o pipefail
    samtools faidx !{ref_cram}
    bazam -n 1  -Xms2G -Xmx20G -Dsamjdk.reference_fasta=!{ref_cram} -bam !{file_tag}.cram | \\
    !{preproc} !{params.bwa_mem} !{ignorealt} !{bwa_opt} -t!{bwa_threads} -R "@RG\\tID:!{file_tag}\\tSM:!{file_tag}\\t!{params.RG}" -p !{ref} - | \\
    !{postalt} samblaster !{samblaster_opt} --addMateTags --ignoreUnmated | \\
    sambamba view -S -f bam -l 0 /dev/stdin | \\
    sambamba sort -t !{sort_threads} -m !{sort_mem}G --tmpdir=!{file_tag}_tmp -o !{file_tag_new}.bam /dev/stdin
    '''
  }
  }
}
if(mode=='fastq'){
  println "fastq mode"
  process fastq_alignment {
    cpus params.cpu
    memory params.mem+'GB'
    tag { "${file_tag}${read_group}" }

  input:
  set val(file_tag), val(nb_groups), val(read_group), file(pair1), file(pair2) from readPairs
  file ref
  file ref_fai
  file ref_sa
  file ref_bwt
  file ref_ann
  file ref_amb
  file ref_pac
  file ref_dict
  file ref_0123
  file ref_bwt8bit
  file ref_alt
	file postaltjs

  output:
	set val(file_tag), val(nb_groups), val(read_group),  file("${file_tag_new}*.bam"), file("${file_tag_new}*.bai") into bam_bai_files0

	//if(!params.recalibration &  !params.input_file){ publishDir "${params.output_folder}/BAM/", mode: 'copy'	}

  shell:
	file_tag_new=file_tag
	bwa_threads  = [params.cpu.intdiv(2) - 1,1].max()
  sort_threads = [params.cpu.intdiv(2) - 1,1].max()
  sort_mem     = [params.mem.intdiv(4),1].max()
  file_tag_new=file_tag_new+"${read_group}"
  if(params.trim) file_tag_new=file_tag_new+'_trimmed'
  if(params.alt)  file_tag_new=file_tag_new+'_alt'
  if(params.alt==null){
          ignorealt='-j'
          postalt=''
  }else{
          ignorealt=''
          postalt="${params.bwakit_root}/k8 ${params.bwakit_root}/bwa-postalt.js "+ref+".alt |"
  }
	if(params.bwa_option_M==null){
          bwa_opt=''
          samblaster_opt=''
        }else{
           bwa_opt='-M '
           samblaster_opt='-M '
        }
  if(nb_groups > 1){
    sort_opt=' -n'
  }else{
    sort_opt=''
  }
	if(params.trim==null){
		'''
    set -o pipefail
    touch !{file_tag_new}.bam.bai
		!{params.bwa_mem} !{ignorealt} !{bwa_opt} -t!{bwa_threads} -R "@RG\\tID:!{file_tag}!{read_group}\\tSM:!{file_tag}\\t!{params.RG}" !{ref} !{pair1} !{pair2} | \\
    !{postalt} samblaster !{samblaster_opt} --addMateTags | \\
    sambamba view -S -f bam -l 0 /dev/stdin | \\
    sambamba sort !{sort_opt} -t !{sort_threads} -m !{sort_mem}G --tmpdir=!{file_tag}_tmp -o !{file_tag_new}.bam /dev/stdin
		'''
 	}else{
		'''
    set -o pipefail
    touch !{file_tag_new}.bam.bai
		AdapterRemoval !{params.adapterremoval_opt} --file1 !{pair1} --file2 !{pair2} --interleaved-output --output1 /dev/stdout | \\
    !{params.bwa_mem} !{ignorealt} !{bwa_opt} -t!{bwa_threads} -R "@RG\\tID:!{file_tag}!{read_group}\\tSM:!{file_tag}\\t!{params.RG}" -p !{ref} - | \\
    !{postalt} samblaster !{samblaster_opt} --addMateTags | \\
    sambamba view -S -f bam -l 0 /dev/stdin | \\
    sambamba sort !{sort_opt} -t !{sort_threads} -m !{sort_mem}G --tmpdir=!{file_tag}_tmp -o !{file_tag_new}.bam /dev/stdin
		'''
	}
     }
}

if(params.input_file){
  bam_bai_files0.into{bam_bai_2group;bam_bai_files2filter}
  bam_bai_grouped4merge = bam_bai_2group.groupTuple(by: 0)
	                                      .map{ a -> [a[0],a[2].size(),a[2],a[3],a[4]] }

	bam_bai_files2filter.filter { a -> a[1] > 1  }
                      .set{mult2QC}

	//QC on each run
	process qualimap_multi {
	    cpus params.cpu
	    memory params.mem+'G'
	    tag { "${file_tag}${read_group}" }

	    publishDir "${params.output_folder}/QC/RG/qualimap/", mode: 'copy'

	    input:
	    set val(file_tag), val(nb_groups), val(read_group),  file(bam), file(bai) from mult2QC
	    file qff from qualimap_ff

	    output:
	    file ("${file_name}") into qualimap_multi_results
	    file ("${file_name}.stats.txt") into flagstat_multi_results

	    shell:
	    feature = qff.name != 'NO_FILE' ? "--feature-file $qff" : ''
      file_name = bam.baseName
	    '''
	    sambamba sort -t !{params.cpu} -m !{params.mem}G --tmpdir=!{file_name}_tmp -o !{file_name}_COsorted.bam !{bam}
	    qualimap bamqc -nt !{params.cpu} !{feature} --skip-duplicated -bam !{file_name}_COsorted.bam --java-mem-size=!{params.mem}G -outdir !{file_name} -outformat html
	    sambamba flagstat -t !{params.cpu} !{bam} > !{file_name}.stats.txt
	    '''
	}

	process multiqc_multi {
	    cpus 2
	    memory '1G'

	    publishDir "${params.output_folder}/QC/RG/qualimap", mode: 'copy'

	    input:
	    file qualimap_results from qualimap_multi_results.collect()
	    file flagstat_results from flagstat_multi_results.collect()
      file multiqc_config from ch_config_for_multiqc

	    output:
	    file("*report.html") into multi_output
	    file("multiqc_multiplex_qualimap_flagstat_report_data/") into multi_output_data

	    shell:
            if( multiqc_config.name=='NO_FILE' ){
                opt = ""
            }else{
                opt = "--config ${multiqc_config}"
            }
	    '''
	    multiqc . -n multiqc_multiplex_qualimap_flagstat_report.html !{opt} --comment "WGS/WES pre-merging QC report"
	    '''
	}

	process merge {
      cpus params.cpu
      memory params.mem+'G'
      tag { file_tag }
      //if(!params.recalibration) publishDir "$params.output_folder/BAM/", mode: 'copy', pattern: "*.bam*"

      input:
      set val(file_tag), val(nb_groups), val(read_group),  file(bams), file(bais) from bam_bai_grouped4merge

      output:
      set val(file_tag), file("${file_tag_new}.bam"), file("${file_tag_new}.bam.bai") into bam_bai_files

      shell:
      file_tag_new=file_tag
      for( rgtmp in read_group ){
        file_tag_new=file_tag_new+"${rgtmp}"
      }
      if(params.trim) file_tag_new=file_tag_new+'_trimmed'
      if(params.alt)  file_tag_new=file_tag_new+'_alt'
      if(nb_groups>1){
         merge_threads  = [params.cpu.intdiv(2) - 1,1].max()
	      sort_threads = [params.cpu.intdiv(2) - 1,1].max()
        sort_mem     = params.mem.div(2)
	      bam_files=" "
	      for( bam in bams ){
        	bam_files=bam_files+" ${bam}"
        }
        file_tag_new=file_tag_new+"_merged"
        if(params.bwa_option_M==null){
	        samblaster_opt=''
	      }else{
	        samblaster_opt='-M '
	      }
        '''
	      sambamba merge -t !{merge_threads} -l 0 /dev/stdout !{bam_files} | \\
        sambamba view -h /dev/stdin | samblaster !{samblaster_opt} --addMateTags | \\
        sambamba view -S -f bam -l 0 /dev/stdin | \\
        sambamba sort -t !{sort_threads} -m !{sort_mem}G --tmpdir=!{file_tag}_tmp -o !{file_tag_new}.bam /dev/stdin
        '''
      }else{
        '''
        touch nomerge
        '''
      }
	}
}

if(params.recalibration){
println "BQSR"

// base quality score recalibration
   process base_quality_score_recalibration {
    cpus params.cpu_BQSR
    memory params.mem_BQSR+'G'
    tag { file_tag }

    //publishDir "$params.output_folder/BAM/", mode: 'copy', pattern: "*bam*"
    publishDir "$params.output_folder/QC/BQSR/", mode: 'copy',
	saveAs: {filename ->
		if (filename.indexOf("table") > 0) "$filename"
		else if (filename.indexOf("plots") > 0) "$filename"
		else null
	}

    input:
    set val(file_tag), file(bam), file(bai) from bam_bai_files
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
    set val(file_tag), file("${file_tag_new}.bam"), file("${file_tag_new}.bam.bai") into final_bam_bai_files, bam2cram

    shell:
    file_name=bam.baseName
    file_tag_new=file_name+'_BQSRecalibrated'
    '''
    gatk BaseRecalibrator --java-options "-Xmx!{params.mem_BQSR}G" -R !{ref} -I !{bam} --known-sites !{known_snps} --known-sites !{known_indels} -O !{file_name}_recal.table
    gatk ApplyBQSR --java-options "-Xmx!{params.mem_BQSR}G" -R !{ref} -I !{bam} --bqsr-recal-file !{file_name}_recal.table -O !{file_tag_new}.bam
    gatk BaseRecalibrator --java-options "-Xmx!{params.mem_BQSR}G" -R !{ref} -I !{file_tag_new}.bam --known-sites !{known_snps} --known-sites !{known_indels} -O !{file_tag_new}_recal.table
    gatk AnalyzeCovariates --java-options "-Xmx!{params.mem_BQSR}G" -before !{file_name}_recal.table -after !{file_tag_new}_recal.table -plots !{file_tag_new}_recalibration_plots.pdf
    mv !{file_tag_new}.bai !{file_tag_new}.bam.bai
    '''
    }
}
else{
  //we duplicate the bam_bai_files into two channels
  bam_bai_files.into {final_bam_bai_files ; bam2cram }
	recal_table_files = Channel.from ( 'NOFILE1', 'NOFILE2' )
}


//These process are always executed

process qualimap_final {
    cpus params.cpu
    memory params.mem+'G'
    tag { file_tag }

    publishDir "${params.output_folder}/QC/qualimap/", mode: 'copy'

    input:
    set val(file_tag), file(bam), file(bai) from final_bam_bai_files
    file qff from qualimap_ff

    output:
    file ("${file_name}") into qualimap_results
    file ("${file_name}.stats.txt") into flagstat_results

    shell:
    feature = qff.name != 'NO_FILE' ? "--feature-file $qff" : ''
    file_name=bam.baseName
    '''
    qualimap bamqc -nt !{params.cpu} !{feature} --skip-duplicated -bam !{bam} --java-mem-size=!{params.mem}G -outdir !{file_name} -outformat html
    sambamba flagstat -t !{params.cpu} !{bam} > !{file_name}.stats.txt
    '''
}

process multiqc_final {
    cpus 2
    memory '2G'

    publishDir "${params.output_folder}/QC/", mode: 'copy'

    input:
    file qualimap_results from qualimap_results.collect()
    file flagstat_results from flagstat_results.collect()
    file BQSR_results from recal_table_files.collect()
    file multiqc_config from ch_config_for_multiqc

    output:
    file("*report.html") into final_output
    file("multiqc_qualimap_flagstat_BQSR_report_data/") into final_output_data

    shell:
    if( multiqc_config.name=='NO_FILE' ){
	opt = ""
    }else{
	opt = "--config ${multiqc_config}"
    }
    '''
    multiqc . -n multiqc_qualimap_flagstat_BQSR_report.html !{opt} --comment "WGS/WES final QC report"
    '''
}


process convert_to_cram {

cpus 4
memory '10G'
def ext = "cram"
def ext_index = "crai"
//we generate the final output
if(params.output_type == "cram") {
  publishDir "${params.output_folder}/CRAM/", mode: 'copy'
}else{
  publishDir "${params.output_folder}/BAM/", mode: 'copy'
  ext = "bam"
  ext_index = "bai"
}


input:
set val(file_tag), file(bam), file(bai) from bam2cram
file(ref) from ch_ref
output:
set val(file_tag), file("${file_name}.${ext}"), file("${file_name}.${ext}.${ext_index}")
script:
file_name=bam.baseName
//BAM->CRAM conversion
if(params.output_type == "cram"){
  """
  samtools view -C  -T ${ref} ${bam} -o ${file_name}.cram
  samtools index ${file_name}.cram
  """
}else{
  """
  """
}


}

//this use ANSI colors to make a short tool description
//useful url: http://www.lihaoyi.com/post/BuildyourownCommandLinewithANSIescapecodes.html
def tool_header (){
        return """
        B\u001b[31;1mW\u001b[33;1mA\u001b[31;1m-MEM\u001b[33;1m : Whole\u001b[32;1m Genome/Exome\u001b[33;1m Alignment\u001b[31;1m (v${workflow.manifest.version})
        """
}


def IARC_Header (){
     return  """
#################################################################################
# ██╗ █████╗ ██████╗  ██████╗██████╗ ██╗ ██████╗ ██╗███╗   ██╗███████╗ ██████╗  #
# ██║██╔══██╗██╔══██╗██╔════╝██╔══██╗██║██╔═══██╗██║████╗  ██║██╔════╝██╔═══██╗ #
# ██║███████║██████╔╝██║     ██████╔╝██║██║   ██║██║██╔██╗ ██║█████╗  ██║   ██║ #
# ██║██╔══██║██╔══██╗██║     ██╔══██╗██║██║   ██║██║██║╚██╗██║██╔══╝  ██║   ██║ #
# ██║██║  ██║██║  ██║╚██████╗██████╔╝██║╚██████╔╝██║██║ ╚████║██║     ╚██████╔╝ #
# ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═════╝ ╚═╝ ╚═════╝ ╚═╝╚═╝  ╚═══╝╚═╝      ╚═════╝  #
# Nextflow pilelines for cancer genomics.########################################
"""
}
