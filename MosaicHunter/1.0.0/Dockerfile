FROM --platform=linux/amd64 ubuntu:latest

RUN apt-get update && \
    apt-get -y install openjdk-8-jdk-headless && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* 

RUN apt-get update
RUN apt-get install -y git
RUN apt-get install -y git-lfs
RUN git lfs install
RUN apt-get install libcurl4

RUN apt-get install -y rsync
RUN rsync -aP hgdownload.soe.ucsc.edu::genome/admin/exe/linux.x86_64/blat/ ./../usr/bin

RUN git clone https://github.com/zzhang526/MosaicHunter
