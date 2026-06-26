# PyPGx benchmark image  ->  pgxbench/pypgx:0.25.0
#   docker build -f benchmark/docker/pypgx.Dockerfile -t pgxbench/pypgx:0.25.0 .
#
# Docs: https://pypgx.readthedocs.io/en/latest/cli.html
#       https://github.com/sbslee/pypgx
# `pypgx run-ngs-pipeline GENE <out_dir> --variants in.vcf.gz [--assembly GRCh37|GRCh38]`
# writes results.zip (SampleTable[Results]) whose data.tsv carries Genotype/Phenotype.
FROM python:3.10-slim

# build toolchain for C-extension pip wheels (pytabix/mappy/pysam)
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential gcc g++ make zlib1g-dev libbz2-dev liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y --no-install-recommends \
        samtools bcftools tabix unzip \
    && rm -rf /var/lib/apt/lists/*

# fiona/pysam wheels pull these; pin pypgx exactly for reproducible allele defs.
RUN pip install --no-cache-dir pypgx==0.25.0 \
    && pypgx --version

# run_pypgx.sh is mounted by Nextflow at runtime; do NOT COPY it here.
