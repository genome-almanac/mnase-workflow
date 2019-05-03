#
# Dockerfile containing v1.10 of cutadapt
# this is the final version to support csqual and csfasta inputs
#

FROM alpine:latest

# scripts
COPY src/ /app

RUN apk add python py-pip python-dev build-base && pip install biopython && \
	apk del py-pip python-dev build-base
