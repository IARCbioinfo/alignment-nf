################## BASE IMAGE #####################
FROM continuumio/miniconda3:4.7.12

################## METADATA #######################

LABEL base_image="continuumio/miniconda3"
LABEL version="4.7.12"
LABEL software="alignment-nf"
LABEL software.version="1.2"
LABEL about.summary="Container image containing all requirements for alignment-nf"
LABEL about.home="http://github.com/IARCbioinfo/alignment-nf"
LABEL about.documentation="http://github.com/IARCbioinfo/alignment-nf/README.md"
LABEL about.license_file="http://github.com/IARCbioinfo/alignment-nf/LICENSE.txt"
LABEL about.license="GNU-3.0"


################## MAINTAINER ######################
MAINTAINER **nalcala** <**alcalan@fellows.iarc.fr**>


################## INSTALLATION ######################
COPY environment.yml /
RUN apt-get update && apt-get install -y procps && apt-get clean -y
RUN conda config --set channel_priority strict
RUN conda env create -n alignment-nf -f /environment.yml && conda clean -a
RUN ln -s /opt/conda/pkgs/bwakit-0.7.15-1/share/bwakit-0.7.15-1/k8 /usr/local/bin/.
RUN wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.0/bwa-mem2-2.0_x64-linux.tar.bz2 && \
    tar jxf bwa-mem2-2.0_x64-linux.tar.bz2 && \
    cp bwa-mem2-2.0_x64-linux/* /usr/local/bin/.
ENV PATH /opt/conda/envs/alignment-nf/bin:$PATH