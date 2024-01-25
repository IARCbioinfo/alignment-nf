################## BASE IMAGE #####################
FROM continuumio/miniconda3:22.11.1

################## METADATA #######################

LABEL base_image="docker pull continuumio/miniconda3"
LABEL version="22.11.1"
LABEL software="alignment-nf"
LABEL software.version="2.1"
LABEL about.summary="Container image containing all requirements for alignment-nf DSL2"
LABEL about.home="http://github.com/IARCbioinfo/alignment-nf"
LABEL about.documentation="http://github.com/IARCbioinfo/alignment-nf/README.md"
LABEL about.license_file="http://github.com/IARCbioinfo/alignment-nf/LICENSE.txt"
LABEL about.license="GNU-3.0"


################## MAINTAINER ######################
LABEL org.opencontainers.image.authors="cahaisv@iarc.who.int"

################## INSTALLATION ######################
COPY environment.yml /
RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps && apt-get clean -y
RUN conda config --set channel_priority strict
RUN conda env create -n alignment-nf -f /environment.yml && conda clean -a

ENV PATH /opt/conda/envs/alignment-nf/bin:$PATH
