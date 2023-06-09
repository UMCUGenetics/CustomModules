FROM --platform=linux/amd64 continuumio/miniconda3

LABEL base_image="continuumio/miniconda3"
LABEL version="1.0.0"
LABEL extra.binaries="vcftools, plink, king"

RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y vcftools=0.1.14 plink=1.90b4 king=2.2.7

