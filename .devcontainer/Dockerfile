# Minimal Docker image using Ubuntu base
FROM ubuntu:22.04
MAINTAINER Niema Moshiri <niemamoshiri@gmail.com>
RUN apt-get update && apt-get -y upgrade && \
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y cmake g++ gcc git make python3 python3-pip wget && \
    python3 -m pip install --upgrade pip && \
    pip3 install bioext && \

    # install cawlign
    git clone "https://github.com/veg/cawlign.git" && \
    cd cawlign && \
    cmake . && \
    make && \
    make install && \
    cd .. && \
    rm -rf cawlign && \

    # Clean up
    rm -rf /root/.cache /tmp/*
