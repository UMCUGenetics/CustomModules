# VCF2GLIMS
Custom vcf parser transforming a (single sample) vcf to GLIMS csv used in Pharmacogenetics workflow.

## Setup development environment
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Build docker container
```bash
docker build . -t umcugenbioinf/vcf2glims:[tag] --no-cache
docker push umcugenbioinf/vcf2glims:[tag]
```

## Build singularity img using dockerhub
```bash
singularity build --disable-cache docker.io-umcugenbioinf-vcf2glims-[tag].img docker://umcugenbioinf/vcf2glims:[tag]
```
