#
# Dockerfile for converting CSFASTA and CSQUAL pairs to FASTQ
#

FROM alpine:latest

# scripts
COPY src/ /app

RUN apk add build-base libstdc++ && cd /app && g++ -std=c++11 convert.cpp -o fastatofastq && \
    mv fastatofastq /bin && cd / && rm -rf /app && apk del build-base
