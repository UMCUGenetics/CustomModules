################## BASE IMAGEs ######################
FROM --platform=linux/amd64 continuumio/miniconda3:24.5.0-0

################## METADATA ######################
LABEL base_image="continuumio/miniconda3:24.5.0-0"
LABEL version="1.1.0"
LABEL extra.binaries="vcftools, plink, king, pandas"

################## INSTALLATIONS ######################
RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y vcftools=0.1.14 plink=1.90b4 king=2.2.7

COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt
