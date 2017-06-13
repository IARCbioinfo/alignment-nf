# Set the base image to Debian
FROM debian:latest

# File Author / Maintainer
MAINTAINER **nalcala** <**alcalan@fellows.iarc.fr**>

RUN apt-get clean && \
  apt-get update -y && \

  # Install dependences
  DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
  make \
  g++ \
  zlib1g-dev \
  libncurses5-dev \
  git \
  wget \
  ca-certificates \
  bzip2 && \

  # Install samtools specific version manually
  wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
  tar -jxf samtools-1.3.1.tar.bz2 && \
  cd samtools-1.3.1 && \
  make && \
  make install && \
  cd .. && \
  rm -rf samtools-1.3.1 samtools-1.3.1.tar.bz2 && \
   
  # Install bwa specific version manually
  wget https://github.com/lh3/bwa/releases/download/v0.7.15/bwakit-0.7.15_x64-linux.tar.bz2 && \
  tar -jxf bwa.kit-0.7.15_x64-linux.tar.bz2 && \
  cp bwa.kit/* /usr/local/bin/.
  rm -rf bwakit-0.7.15_x64-linux.tar.bz2 bwa.kit && \

  # Install samblaster specific version manually
  wget https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.24/samblaster-v.0.1.24.tar.gz && \
  tar -xzf samblaster-v.0.1.24.tar.gz && \
  cd samblaster-v.0.1.24 && \
  make && \
  rm -rf samblaster-v.0.1.24.tar.gz samblaster-v.0.1.24 && \

  # Install sambamba specific version manually
  wget https://github.com/lomereiter/sambamba/releases/download/v0.6.6/sambamba_v0.6.6_linux.tar.bz2 && \
  tar -jxf sambamba_v0.6.6_linux.tar.bz2 && \
  mv sambamba_v0.6.6 /usr/local/bin/sambamba && \
  rm -rf sambamba_v0.6.6_linux.tar.bz2 && \

  # Remove unnecessary dependences
  DEBIAN_FRONTEND=noninteractive apt-get remove -y \
  git && \

  # Clean
  DEBIAN_FRONTEND=noninteractive apt-get autoremove -y && \
  apt-get clean
