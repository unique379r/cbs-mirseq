FROM ubuntu:18.04

## How to built image tis image: docker build . -t cbs-mirseq:v1.0

MAINTAINER Rupesh Kesharwani <bioinforupesh2009.au@gmail.com>

LABEL \
    description="This image is part of CBS-miRSeq workflow."

USER root

ENV TERM xterm-256color

ENV DEBIAN_FRONTEND noninteractive

## gnupg is needed to add new key
RUN apt-get update && apt-get install -y gnupg2

## Install Java 
# Install OpenJDK-8
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk && \
    apt-get install -y ant && \
    apt-get clean;

# Fix certificate issues
RUN apt-get update && \
    apt-get install ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f;

# Install some basic utils
RUN apt-get -y install curl \ 
    unzip \
    gzip \
    gcc \
    g++ \
    git-all \
    make \
    automake \ 
    python-dev \ 
    vim \
    wget \
    tcsh \
    bzip2 \
    libncurses5-dev \
    zlib1g-dev \
    python \
    python-pip

# Install perl and other modules 
RUN apt-get -y install -f cpanminus \
    libpng-dev \
    dos2unix \
    libncurses5-dev \
    libpdf-api2-perl \
    libpdf-api2-simple-perl \
    libreadline6-dev \
    libxt-dev \
    ttf-mscorefonts-installer \
    libbz2-dev \
    libg2-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    r-cran-xml \
    libssl-dev

# Install Bowtie 
ADD https://downloads.sourceforge.net/project/bowtie-bio/bowtie/1.0.0/bowtie-1.0.0-linux-x86_64.zip /tmp
RUN cd /tmp && unzip bowtie-1.0.0-linux-x86_64.zip;
RUN cd /tmp/bowtie-1.0.0 && cp -rav bowtie bowtie-build bowtie-inspect /usr/local/bin

# Install Bedtools
ADD http://github.com/arq5x/bedtools2/releases/download/v2.23.0/bedtools-2.23.0.tar.gz /tmp
RUN cd /tmp && tar -zxvf bedtools-2.23.0.tar.gz
RUN cd /tmp/bedtools2 && make clean && make all && cp -rav bin/* /usr/local/bin

# Install samtools
ADD https://downloads.sourceforge.net/project/samtools/samtools/0.1.18/samtools-0.1.18.tar.bz2 /tmp
RUN cd /tmp && tar -xvjpf samtools-0.1.18.tar.bz2
RUN cd /tmp/samtools-0.1.18 && make && cp -rav samtools /usr/local/bin
RUN cd /tmp/samtools-0.1.18 && cp -rav misc/md5fa misc/md5sum-lite misc/wgsim misc/*.pl bcftools/bcftools /usr/local/bin

# Install samstats
ADD https://downloads.sourceforge.net/project/samstat/samstat-1.5.1.tar.gz /tmp
RUN cd /tmp && tar -zxvf samstat-1.5.1.tar.gz
RUN cd /tmp/samstat-1.5.1 && ./configure && make && make install && cp -v src/samstat /usr/local/bin

# Install fastqc
ADD http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.10.1.zip /tmp
RUN cd /tmp && unzip fastqc_v0.10.1.zip
RUN cd /tmp/FastQC && chmod -R 755 fastqc && cp -av fastqc /usr/local/bin

# Install fastx
ADD http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 /tmp
RUN cd /tmp && tar -xvjpf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 && mv ./bin ./fastx_toolkit_0.0.13 && cp -av fastx_toolkit_0.0.13/* /usr/local/bin

# Install cutadapt 1.3
ADD https://www.dropbox.com/s/l9u22kr8zbgi3ky/cutadapt-1.3.tar.gz?dl=1 /tmp
RUN cd /tmp && tar -zxvf cutadapt-1.3.tar.gz
RUN cd /tmp/cutadapt-1.3 && python setup.py build && python setup.py install && cp -av bin/cutadapt /usr/local/bin

# Install featureCounts v1.4.6 
ADD https://downloads.sourceforge.net/project/subread/subread-1.4.6-p1/subread-1.4.6-p1-Linux-x86_64.tar.gz /tmp
RUN cd /tmp && tar -zxvf subread-1.4.6-p1-Linux-x86_64.tar.gz && cp -av subread-1.4.6-p1-Linux-x86_64/bin/featureCounts /usr/local/bin

# Install RNAhybrid v2.1.1
ADD https://www.dropbox.com/s/qvpb4yu7agjitby/RNAhybrid-2.1.1-src.tar.gz?dl=1 /tmp
RUN cd /tmp/ && tar -xvf RNAhybrid-2.1.1-src.tar.gz
RUN cd /tmp/RNAhybrid-2.1.1 && ./configure && make && make install
RUN cd /tmp/RNAhybrid-2.1.1 && cp -av src/RNAcalibrate src/RNAhybrid src/RNAeffective /usr/local/bin

# Install miranda v3.3a
ADD http://cbio.mskcc.org/microrna_data/miRanda-aug2010.tar.gz /tmp
RUN cd /tmp && tar -zxvf miRanda-aug2010.tar.gz
RUN cd /tmp/miRanda-3.3a && ./configure && make && make install && cp -av src/miranda /usr/local/bin

# Install bfasta/solid2fastq
ADD https://downloads.sourceforge.net/project/bfast/bfast/0.6.5/bfast-0.6.5a.tar.gz /tmp
RUN cd /tmp && tar -zxvf bfast-0.6.5a.tar.gz && cd bfast-0.6.5a && sh autogen.sh && ./configure && make && cp -av scripts/solid2fastq /usr/local/bin

# Install R
RUN apt-get install -y apt-transport-https software-properties-common
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
RUN apt-get update
RUN apt-get install -y r-base
RUN apt-get install -y build-essential

## Copy CBS-miRSeq Package
## It is required to install all R packages before to excute any module of cbs-mirseq
## Just hit: docker run -ti cbs-mirseq:v1.0
## Rscript CBS-miRSeq.v1.0/INSTALL/CBS-miRSeq.Required.PackagesV1.1.0.R
ADD CBS-miRSeq.v1.0 /CBS-miRSeq.v1.0


