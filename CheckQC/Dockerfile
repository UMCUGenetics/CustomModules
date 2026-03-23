################## BASE IMAGE ######################
FROM --platform=linux/amd64 python:3.9

################## METADATA ######################
LABEL base_image="python:3.9"
LABEL version="1.0.0"
LABEL extra.binaries="pandas,PyYAML,pytest"

################## INSTALLATIONS ######################
COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt