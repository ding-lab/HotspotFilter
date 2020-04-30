# Adapted from https://hub.docker.com/r/zlskidmore/varscan/dockerfile
# and https://medium.com/@chadlagore/conda-environments-with-docker-82cdc9d25754
FROM continuumio/miniconda3

MAINTAINER Matthew Wyczalkowski <m.wyczalkowski@wustl.edu>

RUN conda create -n env python=3.6
RUN echo "source activate env" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH

WORKDIR /opt
# https://github.com/geerlingguy/ansible-role-java/issues/64
RUN mkdir -p /usr/share/man/man1    

# run update
RUN apt-get update -y && apt-get install -y \
    curl \
    libnss-sss \
    vim 

RUN conda install -c bioconda bedtools

# https://github.com/ding-lab/BreakPointSurveyor
#RUN git clone --recursive https://github.com/ding-lab/BreakPointSurveyor.git

COPY . /opt/HotspotFilter


