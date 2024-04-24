# Minimal Docker image using Ubuntu base
FROM ubuntu:22.04
MAINTAINER Niema Moshiri <niemamoshiri@gmail.com>
RUN apt-get update && apt-get -y upgrade && \
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y cmake g++ gcc git make wget && \

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
