# Dockerfile to build V2DL3 for Eventdisplay

FROM python:3.10-slim
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get install -y --no-install-recommends wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV PATH /conda/bin:$PATH
ENV PYTHONPATH=$PYTHONPATH:"/V2DL3"

COPY V2DL3 /V2DL3

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /conda && \
    rm ~/miniconda.sh && \
    /conda/bin/conda clean -tipsy && \
    /conda/bin/conda install -c conda-forge mamba && \
    cd /V2DL3 && \
    mamba env update -n base --file environment-eventdisplay.yml && \
    mamba clean --all && conda remove --yes mamba && conda clean --all

SHELL ["/bin/bash", "-c"]

RUN mkdir /data/
RUN source /root/.bashrc && \
    conda init bash
