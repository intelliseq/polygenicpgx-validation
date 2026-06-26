# PGxPOP benchmark image  ->  pgxbench/pgxpop:1.0
#   docker build -f benchmark/docker/pgxpop.Dockerfile -t pgxbench/pgxpop:1.0 .
#
# Docs/repo: https://github.com/PharmGKB/PGxPOP   (a.k.a. Helix-Research-Lab/PGxPOP)
#            allele definitions = PharmCAT (frozen April 2020; project unmaintained).
# PGxPOP calls star alleles + CPIC/DPWG phenotypes from a tabix-indexed VCF, GRCh38 by
# default (it bundles GRCh38 definitions; --build hg19 would liftover). CLI:
#   python bin/PGxPOP.py --vcf in.vcf.gz -g CYP2C19 [--phased] --build grch38 -o out.tsv
#   Output is a TSV with columns:
#     sample_id  gene  diplotype  hap_1  hap_2  hap_1_function  hap_2_function
#     hap_1_variants  hap_2_variants  phenotype  ...  activity_score  uncallable ...
#   Gene set: CFTR CYP2B6 CYP2C9 CYP2C19 CYP2D6 CYP3A5 CYP4F2 DPYD IFNL3 NUDT15
#             SLCO1B1 TPMT UGT1A1 VKORC1
#
# requirements.txt pins: numpy==1.15.2, pyliftover==0.4, pytabix==0.1  -> needs py3.6-3.8.
FROM python:3.8-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        tabix bcftools git build-essential zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Clone the repo (no version tags published -> pin to the default branch HEAD at build).
# numpy 1.15.2 won't build on py3.8; relax to a compatible 1.x and add pyliftover/pytabix.
RUN git clone --depth 1 https://github.com/PharmGKB/PGxPOP.git /opt/PGxPOP \
 && pip install --no-cache-dir "numpy>=1.17,<1.25" "pyliftover==0.4" "pytabix==0.1" \
 && python3 /opt/PGxPOP/bin/PGxPOP.py --help >/dev/null

ENV PGXPOP_DIR=/opt/PGxPOP
# run_pgxpop.sh is mounted by Nextflow at runtime; do NOT COPY it here.
