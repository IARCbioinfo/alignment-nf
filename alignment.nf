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
params.snp_vcf      = "dbsnp.vcf"
params.indel_vcf    = "Mills_1000G_indels.vcf"
params.postaltjs    = "/opt/conda/envs/alignment-nf/share/bwakit-0.7.15-1/bwa-postalt.js"
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

log.info ""
log.info "--------------------------------------------------------"
log.info "  alignment-nf 1.3.0: alignment/realignment workflow for whole exome/whole genome sequencing "
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
    log.info '--input_file     STRING              Input file (comma-separated) with 4 columns:'
    log.info '                                     SM(sample name), RG (read_group_ID), pair1 (path to fastq pair 1), '
    log.info '                                     and pair2 (path to fastq pair 2).'
    log.info '--output_folder  STRING              Output folder (default: .).'
    log.info '--cpu            INTEGER             Number of cpu used by bwa mem and sambamba (default: 8).'
    log.info '--mem            INTEGER             Size of memory used for alignment (in GB) (default: 32).'
    log.info '--RG             STRING              Samtools read group specification with "\t" between fields.'
    log.info '                                         e.g. --RG "PL:ILLUMINA\tDS:custom_read_group".'
    log.info '                                         Default: "PL:ILLUMINA".'
    log.info '--fastq_ext      STRING              Extension of fastq files (default: fastq.gz)'
    log.info '--suffix1        STRING              Suffix of fastq files 1 (default : _1)'
    log.info '--suffix2        STRING              Suffix of fastq files 2 (default : _2)'
    log.info '--snp_vcf        STRING              Path to SNP VCF from GATK bundle (default: dbsnp.vcf)'
    log.info '--indel_vcf      STRING              Path to indel VCF from GATK bundle (default: Mills_1000G_indels.vcf)'
    log.info '--postaltjs      STRING              Path to postalignment javascript bwa-postalt.js'
    log.info '--feature_file   STRING              Path to feature file for qualimap (default: NO_FILE)'
    log.info '--mem_BQSR       INTEGER             Size of memory used for GATK BQSR (in GB) (default: 10)'
    log.info '--cpu_BQSR       INTEGER             Number of cpu used by GATK BQSR (default: 2)'
    log.info '--multiqc_config STRING              Config yaml file for multiqc (default : none)'
    log.info '--adapterremoval_opt STRING          Command line options for AdapterRemoval (default : none)'
    log.info ""
    log.info "Flags:"
    log.info '--trim                               Enable adapter sequence trimming'
    log.info '--recalibration                      Performs base quality score recalibration (GATK)'
    log.info '--alt                                Enable alternative contig handling (for reference genome hg38)'
    log.info '--bwa_option_M                       Trigger the -M option in bwa and the corresponding compatibility option in samblaster'
    log.info '--bwa_mem                            Use "bwa mem" command instead of default: "bwa-mem2 mem"'
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
  log.info "snp_vcf=${params.snp_vcf}"
  log.info "indel_vcf=${params.indel_vcf}"
  log.info "postaltjs=${params.postaltjs}"
  log.info "feature_file=${params.feature_file}"
  log.info "mem_BQSR=${params.mem_BQSR}"
  log.info "cpu_BQSR=${params.cpu_BQSR}"
  log.info "multiqc_config=${params.multiqc_config}"
  log.info "bwa_mem=${params.bwa_mem}"
  log.info "adapterremoval_opt=${params.adapterremoval_opt}"
  log.info "recalibration=${params.recalibration}"
  log.info "alt=${params.alt}"
  log.info "trim=${params.trim}"
  log.info "bwa_option_M=${params.bwa_option_M}"
  log.info "help=${params.help}"
}






/***************************************************************************************/
/************************ handle global parameters *************************************/
/***************************************************************************************/

ignorealt = params.alt ? '': '-j'
postalt   = params.alt ? 'k8 bwa-postalt.js '+ file(params.ref+'.alt').name + ' | ' : ''
bwa_opt   = params.bwa_option_M ? '-M ' : ''
samblaster_opt= params.bwa_option_M ? '-M ' : ''

//multiqc config file
ch_config_for_multiqc = file(params.multiqc_config)

//reference file and its indexes for bwa
bwa_ref = tuple file(params.ref), file(params.ref+'.fai'), 
  file(params.ref+'.sa'), file(params.ref+'.bwt'), 
  file(params.ref+'.ann'), file(params.ref+'.amb'), file(params.ref+'.pac'), 
  file(params.ref.replaceFirst(/fasta/, "").replaceFirst(/fa/, "") +'dict'), 
  params.bwa_mem != "bwa-mem2 mem" ? file('NO_0123') : file(params.ref+'.0123'), 
  params.bwa_mem != "bwa-mem2 mem" ? file('NO_bwtnbit') : file(params.ref+'*bwt.*bit.*'),
  params.alt ? file(params.ref+'.alt') : file('NO_ALT')

println("${bwa_ref[9]}")

postaltjs = file( params.postaltjs )

//reference file and its indexes for baserecalibration
bqsr_ref = tuple file(params.ref), file(params.ref+'.fai'), file(params.ref.replaceFirst(/fasta/, "").replaceFirst(/fa/, "") +'dict')

//get know site VCFs from GATK bundle
known_snps         = tuple file( params.snp_vcf ), file( params.snp_vcf+'.tbi' )
known_indels       = tuple file( params.indel_vcf ), file( params.indel_vcf+'.tbi' )

//qualimap feature file
qualimap_ff = file(params.feature_file)





/***************************************************************************************/
/************************  Process : fastq_alignment ***********************************/
/***************************************************************************************/


process fastq_alignment {

  cpus params.cpu
  memory params.mem+'GB'    

  if(!params.recalibration &  !params.input_file){ publishDir "${params.output_folder}/BAM/", mode: 'copy'	}

  input:
    tuple val(file_tag), val(nb_groups), val(read_group), path(pair1), path(pair2)
    tuple path(ref), path(ref_fai), path(ref_sa), path(ref_bwt), path(ref_ann), path(ref_amb), path(ref_pac), path(ref_dict), path(ref_0123), path(ref_bwt8bit), path(ref_alt)
    path(postaltjs)
                 
  output:
	  tuple val(file_tag), val(nb_groups), val(read_group), path("${file_tag_new}*.bam"), path("${file_tag_new}*.bai"), emit: bamfiles

  shell:
    bwa_threads  = [params.cpu.intdiv(2) - 1,1].max()
    sort_threads = [params.cpu.intdiv(2) - 1,1].max()
    sort_mem     = [params.mem.intdiv(4),1].max()
    file_tag_new = read_group == "" ? file_tag : "${file_tag}_${read_group}"
    if(params.trim) file_tag_new=file_tag_new+'_trimmed'
    if(params.alt)  file_tag_new=file_tag_new+'_alt'	
    sort_opt = nb_groups > 1 ? ' -n' : ''
    RG = "\"@RG\\tID:$file_tag_new\\tSM:$file_tag\\t$params.RG\""

    if(params.trim==null){
      """
      set -o pipefail
      touch ${file_tag_new}.bam.bai
      $params.bwa_mem $ignorealt $bwa_opt -t$bwa_threads -R $RG $ref $pair1 $pair2 | \
      $postalt samblaster $samblaster_opt --addMateTags | \
      sambamba view -S -f bam -l 0 /dev/stdin | \
      sambamba sort $sort_opt -t $sort_threads -m ${sort_mem}G --tmpdir=${file_tag}_tmp -o ${file_tag_new}.bam /dev/stdin
      """
    }else{
      """
      set -o pipefail
      touch ${file_tag_new}.bam.bai
      AdapterRemoval $params.adapterremoval_opt --file1 $pair1 --file2 $pair2 --interleaved-output --output1 /dev/stdout | \
      $params.bwa_mem $ignorealt $bwa_opt -t$bwa_threads -R $RG -p ${ref} - | \
      $postalt samblaster $samblaster_opt --addMateTags | \
      sambamba view -S -f bam -l 0 /dev/stdin | \
      sambamba sort $sort_opt -t $sort_threads -m ${sort_mem}G --tmpdir=${file_tag}_tmp -o ${file_tag_new}.bam /dev/stdin
      """
    }

}

/***************************************************************************************/
/************************  Process : bam_realignment ***********************************/
/***************************************************************************************/

process bam_realignment {
  cpus params.cpu
  memory params.mem+'G'
        
  if(!params.recalibration) publishDir "${params.output_folder}/BAM/", mode: 'copy'

	input:
    path infile
    tuple path(ref), path(ref_fai), path(ref_sa), path(ref_bwt), path(ref_ann), path(ref_amb), path(ref_pac), path(ref_dict), path(ref_0123), path(ref_bwt8bit), path(ref_alt)
    path postaltjs
     
  output:
	  tuple val(file_tag), file("${file_tag_new}*.bam"), file("${file_tag_new}*.bai")

  shell:
	  file_tag = infile.baseName
	  file_tag_new=file_tag+'_realigned'
	  if(params.trim) file_tag_new=file_tag_new+'_trimmed'
	  if(params.alt)  file_tag_new=file_tag_new+'_alt'
    preproc = params.trim ? "AdapterRemoval $params.adapterremoval_opt --interleaved --file1 /dev/stdin --output1 /dev/stdout |" : ''
    bwa_threads  = [params.cpu.intdiv(2) - 1,1].max()
    sort_threads = [params.cpu.intdiv(2) - 1,1].max()
    sort_mem     = params.mem.div(4)
    read_group = "@RG\tID:$file_tag$read_group\tSM:$file_tag\t$params.RG"
    """
    set -o pipefail
    samtools collate -uOn 128 ${file_tag}.bam tmp_$file_tag | \
    samtools fastq - | \
    $preproc $params.bwa_mem $ignorealt $bwa_opt -t$bwa_threads -R $read_group -p $ref - | \
    $postalt samblaster $samblaster_opt --addMateTags --ignoreUnmated | \
    sambamba view -S -f bam -l 0 /dev/stdin | \
    sambamba sort -t $sort_threads -m ${sort_mem}G --tmpdir=${file_tag}_tmp -o ${file_tag_new}.bam /dev/stdin
    """
}

/***************************************************************************************/
/************************  Process : qualimap_multi ************************************/
/***************************************************************************************/

process qualimap_multi {
	
  cpus params.cpu
	memory params.mem+'G'

	publishDir "${params.output_folder}/QC/BAM/qualimap/", mode: 'copy'

	input:
	  tuple val(file_tag), val(nb_groups), val(read_group), path(bam), path(bai)
	  path(qff)

	output:
	  tuple path("${file_name}"), path("${file_name}.stats.txt")

	shell:
	  feature = qff.name != 'NO_FILE' ? "--feature-file $qff" : ''
    file_name = bam.baseName
	  """
	  sambamba sort -t $params.cpu -m ${params.mem}G --tmpdir=${file_name}_tmp -o ${file_name}_COsorted.bam $bam
	  qualimap bamqc -nt $params.cpu $feature --skip-duplicated -bam ${file_name}_COsorted.bam --java-mem-size=${params.mem}G -outdir $file_name -outformat html
	  sambamba flagstat -t $params.cpu $bam > ${file_name}.stats.txt
	  """
}

/***************************************************************************************/
/************************  Process : multiqc_multi *************************************/
/***************************************************************************************/

process multiqc_multi {
	
  cpus 2
	memory '1G'

	publishDir "${params.output_folder}/QC/BAM/qualimap/", mode: 'copy'

	input:
	  path qualimap_results
    path multiqc_config

	output:
    path("*report.html")
	  path("multiqc_*_data/*")

	shell:
    opt = (multiqc_config.name=='NO_FILE' ) ?  "" : "--config ${multiqc_config}"
	  """
	  multiqc . -n multiqc_multiplex_qualimap_flagstat_report.html $opt --comment "WGS/WES pre-merging QC report"
	  """
}


/***************************************************************************************/
/************************  Process : merge *********************************************/
/***************************************************************************************/

process merge_bam {
      
  cpus params.cpu
  memory params.mem+'G'

  if(params.recalibration) publishDir "$params.output_folder/BAM/", mode: 'copy', pattern: "*.bam*"

  input:
    tuple val(file_tag), val(nb_groups), val(read_group), path(bams), path(bais)

  output:
    tuple val(file_tag), path("${file_tag_new}.bam"), path("${file_tag_new}.bam.bai")

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
	    /*for( bam in bams ){
        bam_files=bam_files+" ${bam}"
      }*/
      file_tag_new=file_tag_new+"_merged"
      """
	    sambamba merge -t $merge_threads -l 0 /dev/stdout $bam_files | \
      sambamba view -h /dev/stdin | \
      samblaster $samblaster_opt --addMateTags | \
      sambamba view -S -f bam -l 0 /dev/stdin | \
      sambamba sort -t $sort_threads -m ${sort_mem}G --tmpdir=${file_tag}_tmp -o ${file_tag_new}.bam /dev/stdin
      """
    }else{
      """
      touch nomerge
      """
    }
}

/***************************************************************************************/
/************************  Process : mebase_quality_score_recalibrationrge *************/
/***************************************************************************************/

process base_quality_score_recalibration {

  cpus params.cpu_BQSR
  memory params.mem_BQSR+'G'
    
  publishDir "$params.output_folder/BAM/", mode: 'copy', pattern: "*bam*"
  publishDir "$params.output_folder/QC/BAM/BQSR/", mode: 'copy',
	saveAs: {filename -> 
		if (filename.indexOf("table") > 0) "$filename"
		else if (filename.indexOf("plots") > 0) "$filename"
		else null
	}

  input:
    tuple val(file_tag), path(bam), path(bai)
    tuple path(ref), path(ref_fai), path(ref_dict)
    tuple path(known_snps), path(known_snps_index)
    tuple path(known_indels), path(known_indels_index)
    
  output:
    tuple val(file_tag), path("${file_tag_new}.bam"), path("${file_tag_new}.bam.bai"), emit: bamfiles
    path("*_recal.table"), emit: recal_table_files
    path("*plots.pdf")

  shell:
    file_name=bam.baseName
    file_tag_new=file_name+'_BQSRecalibrated'
    """
    gatk BaseRecalibrator --java-options "-Xmx${params.mem_BQSR}G" -R $ref -I $bam --known-sites $known_snps --known-sites $known_indels -O ${file_name}_recal.table
    gatk ApplyBQSR --java-options "-Xmx${params.mem_BQSR}G" -R $ref -I $bam --bqsr-recal-file ${file_name}_recal.table -O ${file_tag_new}.bam
    gatk BaseRecalibrator --java-options "-Xmx${params.mem_BQSR}G" -R $ref -I ${file_tag_new}.bam --known-sites ${known_snps} --known-sites ${known_indels} -O ${file_tag_new}_recal.table
    gatk AnalyzeCovariates --java-options "-Xmx${params.mem_BQSR}G" -before ${file_name}_recal.table -after ${file_tag_new}_recal.table -plots ${file_tag_new}_recalibration_plots.pdf
    touch ${file_tag_new}_recalibration_plots.pdf	
    mv ${file_tag_new}.bai ${file_tag_new}.bam.bai
    """
}

/***************************************************************************************/
/************************  Process : qualimap_final ************************************/
/***************************************************************************************/

process qualimap_final {
  
  cpus params.cpu
  memory params.mem+'G'

  publishDir "${params.output_folder}/QC/BAM/qualimap/", mode: 'copy'

  input:
    tuple val(file_tag), path(bam), path(bai)
    path(qff)

  output:
    tuple path("$file_name"), path("${file_name}.stats.txt")

  shell:
    feature = qff.name != 'NO_FILE' ? "--feature-file $qff" : ''
    file_name=bam.baseName
    """
    qualimap bamqc -nt $params.cpu $feature --skip-duplicated -bam $bam --java-mem-size=${params.mem}G -outdir $file_name -outformat html
    sambamba flagstat -t $params.cpu $bam > ${file_name}.stats.txt
    """
}

/***************************************************************************************/
/************************  Process : multiqc_final *************************************/
/***************************************************************************************/

process multiqc_final {
    
  cpus 2
  memory '2G'

  publishDir "${params.output_folder}/QC/BAM/", mode: 'copy'

  input:
    path(qualimap_results)
    path(BQSR_results)
    path(multiqc_config)

  output:
    file("*report.html")
    file("multiqc_*_data/*")

  shell:
    opt = (multiqc_config.name=='NO_FILE') ? "" : "--config ${multiqc_config}"
    """
    multiqc . -n multiqc_qualimap_flagstat_BQSR_report.html $opt --comment "WGS/WES final QC report"
    """
}




/***************************************************************************************/
/************************  Workflow : main *********************************************/
/***************************************************************************************/



workflow {

  mode='fastq'

  if(params.input_file){

    readPairs = Channel.fromPath(params.input_file).splitCsv(header: true, sep: '\t', strip: true) | map {
        row ->
          assert (row.SM != null ) : "Error: SM column is missing, check your input file"
          assert (row.RG != null ) : "Error: RG column is missing, check your input file"
          assert (row.pair1 != null ) : "Error: pair1 column is missing, check your input file"
          assert (row.pair2 != null ) : "Error: pair2 column is missing, check your input file"
          tuple( row.SM, row.RG, file(row.pair1), file(row.pair2) )} | groupTuple(by:0) | map{
            row -> tuple(row[0], row[1].size(), row[1], row[2], row[3])} | transpose()

  }else if(params.input_folder){
    
    if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.fastq_ext}/ }.size() > 0){
      
      println "fastq files found, proceed with alignment"
      readPairs = Channel.fromFilePairs(params.input_folder +"/*{${params.suffix1},${params.suffix2}}" +'.'+ params.fastq_ext)
        .map { row -> tuple( row[0] , 1 , "NO_RG" , row[1][0], row[1][1] )}

    }else if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0){
      
      println "BAM files found, proceed with realignment"
      files = Channel.fromPath( params.input_folder+'/*.bam' )
      mode='bam'
    
    }else{

      println "ERROR: input folder contains no fastq nor BAM files"; System.exit(0)

    }
  }

  // bwa alignment
  if(mode=="fastq"){
    fastq_alignment(readPairs,bwa_ref,postaltjs)
    bamfiles = fastq_alignment.out.bamfiles

    // separate standalone bam and bam that need to be merged
    bamBranch = bamfiles.groupTuple(by:0)
      .map{ row -> tuple(row[0], row[2].size(), row[2], row[3], row[4])}
      .branch{ single: it[1]==1
               multi: it[1]> 1
      }
    //bamBranch.single.view()
    //bamBranch.multi | view()

    // QC and merge bam
    qualimap_multi( bamfiles.filter{ it[1]>1 }, qualimap_ff)
    multiqc_multi(qualimap_multi.out.collect(), ch_config_for_multiqc)
    multi=merge_bam(bamBranch.multi)

    // gather standalone bam files and merged bam files
    single = bamBranch.single.map{ row -> tuple(row[0], row[3], row[4]) }
    bamfiles = single.concat(multi)

  } else if(mode=="bam"){
    bam_realignment(files,bwa_ref,postaltjs)
    bamfiles = bam_realignment.out.bamfiles
  }
  
  
  // BQSR recalibration
  if(params.recalibration){
    
    base_quality_score_recalibration(bamfiles,bqsr_ref,known_snps,known_indels)
    bamfiles = base_quality_score_recalibration.out.bamfiles
    recal_table_files = base_quality_score_recalibration.out.recal_table_files
  } else {
    recal_table_files = file("NO_BQSR")
  }

  // Qualimap 
  qualimap_final(bamfiles,qualimap_ff)

  // multiqc
  multiqc_final(qualimap_final.out.collect(), recal_table_files, ch_config_for_multiqc)


}




