################## BASE IMAGE ######################
FROM --platform=linux/amd64 python:3.9

################## METADATA ######################
LABEL base_image="python:3.9"
LABEL version="1.0.1"
LABEL extra.binaries="pandas,PyYAML,pytest"

################## INSTALLATIONS ######################
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
RUN python3 -m venv $POETRY_VENV \
	&& $POETRY_VENV/bin/pip3 install -U pip3 setuptools \
	&& $POETRY_VENV/bin/pip3 install poetry==${POETRY_VERSION}

# Copy project requirement files here to ensure they will be cached.
COPY pyproject.toml poetry.lock ./

# Install dependencies.
RUN --mount=type=cache,target=$POETRY_CACHE_DIR poetry install --without dev --no-interaction --no-cache
