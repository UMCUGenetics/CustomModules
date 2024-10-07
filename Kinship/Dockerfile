################## BASE IMAGEs ######################
FROM continuumio/miniconda3:24.5.0-0

################## METADATA ######################
LABEL base_image="continuumio/miniconda3:24.5.0-0"
LABEL version="1.1.1"
LABEL extra.binaries="vcftools, plink, king, pandas"

################## INSTALLATIONS ######################
# Use conda to install other dependencies.
RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y vcftools=0.1.14 plink=1.90b4 king=2.2.7

# Use poetry to install virtualenv.
RUN pip install poetry==1.8.3
ENV POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_IN_PROJECT=1 \
    POETRY_VIRTUALENVS_CREATE=1

COPY pyproject.toml poetry.lock ./

RUN poetry install --without dev --no-interaction --no-cache
ENV VIRTUAL_ENV=/.venv \
    PATH="/.venv/bin:$PATH"
