################## BASE IMAGE #####################
FROM continuumio/miniconda3:4.7.12

################## METADATA #######################

ARG bwa_mem2_version=2.2.1


LABEL base_image="continuumio/miniconda3"
LABEL version="4.7.12"
LABEL software="alignment-nf"
LABEL software.version="1.3"
LABEL about.summary="Container image containing all requirements for alignment-nf"
LABEL about.home="http://github.com/IARCbioinfo/alignment-nf"
LABEL about.documentation="http://github.com/IARCbioinfo/alignment-nf/README.md"
LABEL about.license_file="http://github.com/IARCbioinfo/alignment-nf/LICENSE.txt"
LABEL about.license="GNU-3.0"


################## MAINTAINER ######################
MAINTAINER **nalcala** <**alcalan@fellows.iarc.fr**>


################## INSTALLATION ######################
COPY environment.yml /
RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps && apt-get clean -y
RUN conda config --set channel_priority strict
RUN conda env create -n alignment-nf -f /environment.yml && conda clean -a
RUN ln -s /opt/conda/pkgs/bwakit-0.7.15-1/share/bwakit-0.7.15-1/k8 /usr/local/bin/.
RUN wget -q https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-${bwa_mem2_version}_x64-linux.tar.bz2 && \
    tar jxf bwa-mem2-${bwa_mem2_version}_x64-linux.tar.bz2 && \
    mv bwa-mem2-${bwa_mem2_version}_x64-linux/* /usr/local/bin/. && \
    rm -rf bwa-mem2-${bwa_mem2_version}
ENV PATH /opt/conda/envs/alignment-nf/bin:$PATH
