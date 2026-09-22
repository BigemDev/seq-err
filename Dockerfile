FROM python:3.11-slim

RUN apt-get update -qq && apt-get install -y -qq \
    build-essential zlib1g-dev libbz2-dev liblzma-dev \
    libcurl4-openssl-dev default-jre wget unzip git \
    samtools bcftools tabix bwa >/dev/null

RUN wget -q https://github.com/broadinstitute/gatk/releases/download/4.7.0.0/gatk-4.7.0.0.zip && \
    unzip -q gatk-4.7.0.0.zip && \
    mv gatk-4.7.0.0/gatk /usr/local/bin/ && \
    mv gatk-4.7.0.0/gatk-package-4.7.0.0-local.jar /usr/local/bin/ && \
    rm -rf gatk-4.7.0.0*

RUN pip install --break-system-packages -q cffi && \
    git clone --quiet --recursive https://github.com/nanoporetech/bwapy.git /tmp/bwapy && \
    cd /tmp/bwapy && \
    make bwa/libbwa.a && \
    pip install --break-system-packages -q . && \
    rm -rf /tmp/bwapy

WORKDIR /app
COPY requirements.txt requirements-dev.txt ./
RUN pip install -r requirements-dev.txt --break-system-packages -q