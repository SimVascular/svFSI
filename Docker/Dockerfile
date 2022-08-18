FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

# Install supporting packages
RUN apt-get update
RUN apt-get install -yq --no-install-recommends apt-utils
RUN apt-get install -yq git make cmake
RUN apt-get install -yq libblas-dev liblapack-dev
RUN apt-get install -yq gcc g++ gfortran 
RUN apt-get install -yq openmpi-bin libopenmpi-dev

# Create directory
RUN mkdir -p /home/test

# Add non-root user and set up home directory
RUN useradd testuser -u 1000 -g 100 -m -s /bin/bash
RUN chown testuser /home/test
USER testuser
WORKDIR /home/test

# Obtain source code
RUN git clone https://github.com/SimVascular/svFSI

# Compile svFSI source code
RUN mkdir Build && \
    cd Build && \
    cmake ../svFSI && \
    make

ENV PATH=$PATH:/home/test/Build/svFSI-build/bin
