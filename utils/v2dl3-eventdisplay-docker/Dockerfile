# Dockerfile to build V2DL3 for Eventdisplay
FROM mambaorg/micromamba:latest

RUN micromamba install -y -n base -c conda-forge git python="3.11"

RUN eval "$(micromamba shell hook --shell bash)" && \
    micromamba activate && \
    git clone https://github.com/VERITAS-Observatory/V2DL3.git V2DL3/ && \
    cd V2DL3 && \
    git checkout pypi-eventdisplay

RUN cd V2DL3 && \
    micromamba install -y -n base -f environment-eventdisplay.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1


RUN python -m pip install --upgrade pip

RUN cd V2DL3 && \
    pip install --no-cache-dir .

SHELL ["/bin/bash", "-c"]

RUN mkdir data/
