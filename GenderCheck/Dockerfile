################## BASE IMAGE ######################
FROM --platform=linux/amd64 python:3.11

################## METADATA ######################
LABEL base_image="python:3.11"
LABEL version="1.0.0"
LABEL extra.binaries="pysam,pytest"

################## INSTALLATIONS ######################
COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt