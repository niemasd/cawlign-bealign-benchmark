# Minimal Docker image using Alpine base
FROM alpine:3.13.5
MAINTAINER Niema Moshiri <niemamoshiri@gmail.com>
RUN apk update && \
    apk add bash cmake gcc g++ make && \
    wget -qO- "https://github.com/veg/tn93/archive/refs/tags/v.1.0.11.tar.gz" | tar -zx && \
    cd tn93-* && \
    cmake . && \
    make install && \
    cd .. && \
    rm -rf tn93-*
