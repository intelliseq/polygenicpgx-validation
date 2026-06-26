# Aldy benchmark image  ->  pgxbench/aldy:4.6
#   docker build -f benchmark/docker/aldy.Dockerfile -t pgxbench/aldy:4.6 .
#
# Docs: https://github.com/0xTCG/aldy   https://aldy.readthedocs.io/
# `aldy genotype -p illumina -g GENE --genome hg38|hg19 [-r ref.fa] -o out.aldy in.cram`
# Output .aldy file: a `#Solution N: ...` comment header plus the diplotype in the
# `Major` column using "*1/*4" slash notation. Aldy ships the CBC ILP solver wheel.
FROM python:3.10-slim

# build toolchain for C-extension pip wheels (pytabix/mappy/pysam)
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential gcc g++ make zlib1g-dev libbz2-dev liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y --no-install-recommends \
        samtools bcftools tabix \
    && rm -rf /var/lib/apt/lists/*

# Aldy 4.6 (4.x line); bundled CBC solver is sufficient (no Gurobi license needed).
RUN pip install --no-cache-dir "aldy==4.6" \
    && aldy --help >/dev/null

# run_aldy.sh is mounted by Nextflow at runtime; do NOT COPY it here.
