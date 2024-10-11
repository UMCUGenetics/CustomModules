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
ENV POETRY_VERSION=1.8.3
ENV POETRY_HOME=/opt/poetry
ENV POETRY_VENV=/opt/poetry-venv

# Tell Poetry where to place its cache and virtual environment
ENV POETRY_CACHE_DIR=/tmp/poetry_cache

# Do not ask any interactive question
ENV POETRY_NO_INTERACTION=1

# Make poetry create the virtual environment in the project's root
# it gets named `.venv`
ENV POETRY_VIRTUALENVS_IN_PROJECT=1
ENV POETRY_VIRTUALENVS_CREATE=1

# Set virtual_env variable
ENV VIRTUAL_ENV=/.venv
# Prepend virtual environments path
ENV PATH="${VIRTUAL_ENV}/bin:${POETRY_VENV}/bin:${PATH}"

# Creating a virtual environment just for poetry and install it with pip
RUN ls -l `which pip`
RUN ls -l `which pip3`
RUN python3 -m venv $POETRY_VENV \
	&& $POETRY_VENV/bin/pip install -U pip setuptools \
	&& $POETRY_VENV/bin/pip install poetry==${POETRY_VERSION}

# Copy project requirement files here to ensure they will be cached.
COPY pyproject.toml poetry.lock ./

# Install dependencies.
RUN --mount=type=cache,target=$POETRY_CACHE_DIR poetry install --without dev --no-interaction --no-cache
