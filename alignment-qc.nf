


//we display the header of the tool
log.info IARC_Header()
log.info tool_header()

//we init some basic parameters
params.input_folder = null
params.cpu          = 8
params.mem          = 32
params.help         = null
params.feature_file = 'NO_FILE'
params.multiqc_config = 'NO_FILE'

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ''
    log.info 'nextflow run iarcbioinfo/alignment-qc.nf [-with-docker] --input_folder input/ [OPTIONS]'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '--input_folder   FOLDER              Folder containing BAM or fastq files to be aligned.'
    log.info 'Optional arguments:'
    log.info '--output_folder  STRING              Output folder (default: "./results").'
    log.info '--multiqc_config STRING              Config yaml file for multiqc (default : none)'
    log.info '--feature_file   STRING              Path to feature file for qualimap (default: NO_FILE)'

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

assert (params.input_folder != null) : "please specify input_folder"
mode=''
//we try BAM files
if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0){
    println "BAM files found, proceed with realignment";
    mode ='bam'
    bams = Channel.fromPath( params.input_folder+'/*.bam' )
            .map {path -> [ path.name.replace(".bam",""),path]}
    bams_index = Channel.fromPath( params.input_folder+'/*.bam.bai')
            .map {  path -> [ path.name.replace(".bam.bai",""), path ] }
    //we create the chanel
    final_bam_bai_files = bams.join(bams_index)

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
          final_bam_bai_files = crams.join(crams_index)
    }else{
      	println "ERROR: input folder contains no BAM/CRAM files"; System.exit(1)
    }
  }
//we dont use this channel at the moment
recal_table_files = Channel.from ( 'NOFILE1', 'NOFILE2' )

//multiqc config file
ch_config_for_multiqc = file(params.multiqc_config)
//qualimap feature file
qualimap_ff = file(params.feature_file)

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
  if(mode == "bam"){
    '''
    qualimap bamqc -nt !{params.cpu} !{feature} --skip-duplicated -bam !{bam} --java-mem-size=!{params.mem}G -outdir !{file_name} -outformat html
    samtools flagstat !{bam} > !{file_name}.stats.txt
    '''
   }else{
     //we have to process cram file
     '''
     mkfifo inbam
     samtools view -b !{bam} > inbam &
     qualimap bamqc -nt !{params.cpu} !{feature} --skip-duplicated -bam inbam --java-mem-size=!{params.mem}G -outdir !{file_name} -outformat html
     samtools flagstat !{bam} > !{file_name}.stats.txt
     '''
   }
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








//this use ANSI colors to make a short tool description
//useful url: http://www.lihaoyi.com/post/BuildyourownCommandLinewithANSIescapecodes.html
def tool_header (){
        return """
        Q\u001b[31;1mC\u001b[33;1m-ALN\u001b[33;1m : QC\u001b[32;1m for Genome/Exome\u001b[33;1m Alignments\u001b[31;1m (v${workflow.manifest.version})
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
# Nextflow pipelines for cancer genomics.########################################
"""
}
