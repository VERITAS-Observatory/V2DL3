# Dockerfile to build V2DL3 for Eventdisplay
FROM mambaorg/micromamba:latest
WORKDIR /workdir/

RUN micromamba install -y -n base -c conda-forge git

RUN eval "$(micromamba shell hook --shell bash)" && \
    micromamba activate && \
    git clone --depth 1 https://github.com/VERITAS-Observatory/V2DL3.git V2DL3/

RUN cd V2DL3 && \
    micromamba install -y -n base -f environment-eventdisplay.yml && \
    micromamba clean --all --yes
ENV PYTHONPATH="/workdir/V2DL3/"

SHELL ["/bin/bash", "-c"]

RUN mkdir data/
