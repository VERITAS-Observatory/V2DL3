# Dockerfile to build V2DL3 for Eventdisplay
FROM python:3.11-slim
WORKDIR /workdir/

RUN apt-get update && apt-get install -y \
    python3-pip \
    git \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/VERITAS-Observatory/V2DL3.git V2DL3

RUN cd V2DL3 && \
    pip install --no-cache-dir .

RUN mkdir data/
