# Minimal Docker image using Ubuntu base
FROM ubuntu:22.04
MAINTAINER Niema Moshiri <niemamoshiri@gmail.com>
RUN apt-get update && apt-get -y upgrade && \
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y cmake g++ gcc make wget && \

    # install tn93
    wget -qO- "https://github.com/veg/tn93/archive/refs/tags/v.1.0.13.actual.tar.gz" | tar -zx && \
    cd tn93-* && \
    cmake . && \
    make install && \
    cd .. && \
    rm -rf tn93-*
