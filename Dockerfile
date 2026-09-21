FROM python:3.11-slim

RUN apt-get update -qq && apt-get install -y -qq \
    build-essential zlib1g-dev libbz2-dev liblzma-dev \
    libcurl4-openssl-dev default-jre wget unzip \
    samtools bcftools tabix >/dev/null

RUN wget -q https://github.com/broadinstitute/gatk/releases/download/4.7.0.0/gatk-4.7.0.0.zip && \
    unzip -q gatk-4.7.0.0.zip && \
    mv gatk-4.7.0.0/gatk /usr/local/bin/ && \
    mv gatk-4.7.0.0/gatk-package-4.7.0.0-local.jar /usr/local/bin/ && \
    rm -rf gatk-4.7.0.0*

WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt --break-system-packages -q