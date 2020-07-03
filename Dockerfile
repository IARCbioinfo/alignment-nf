################## BASE IMAGE #####################
FROM continuumio/miniconda3:4.7.12

################## METADATA #######################

LABEL base_image="continuumio/miniconda3"
LABEL version="4.7.12"
LABEL software="alignment-nf"
LABEL software.version="1.1"
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
RUN conda env create -n alignment-nf -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/alignment-nf/bin:$PATH
