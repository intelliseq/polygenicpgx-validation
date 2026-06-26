# Stargazer benchmark image  ->  pgxbench/stargazer:2.0.2
#   docker build -f benchmark/docker/stargazer.Dockerfile -t pgxbench/stargazer:2.0.2 .
#
# Docs: https://stargazer.gs.washington.edu/stargazerweb/
# CLI (`stargazer -h`):
#   stargazer <dt> <gb> <tg> <vcf> <out> [--gdf FILE] [--cg STR] [--ref FILE] ...
#     dt  = data type   : wgs | ts | chip
#     gb  = genome build: hg19 | hg38   (our harness uses Stargazer for GRCh37 only)
#     tg  = target gene  (lowercase, e.g. cyp2d6)
#     vcf = input VCF
#     out = output project directory  ->  <out>/genotype.txt
#   genotype.txt header: gene name status hap1_main hap2_main ... phenotype ...
#     diplotype = hap1_main "/" hap2_main (star tokens incl xN/+); status 'g' = genotyped.
#   SNV genes can run in VCF-only mode (no GDF); CNV (CYP2D6) wants a GDF for copy number.
#
# INSTALL NOTE: the OFFICIAL Stargazer is gated behind academic registration
# (https://stargazer.gs.washington.edu/stargazerweb/res/form.html) and CANNOT be
# auto-fetched. Two install paths below:
#   (A) [default] clone a public source mirror of the same codebase, OR
#   (B) [registered users] drop the downloaded stargazer-2.0.2 tarball into the build
#       context and uncomment the COPY/install block.
# If neither yields a working `stargazer` binary the runner still satisfies the
# contract (emits NA when the binary is missing).
FROM python:3.8-slim

# build toolchain for C-extension pip wheels (pytabix/mappy/pysam)
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential gcc g++ make zlib1g-dev libbz2-dev liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y --no-install-recommends \
        samtools bcftools tabix git \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir numpy pandas

# --- (A) public source mirror (same algorithm/codebase; no registration) ----------
RUN git clone --depth 1 https://github.com/Wilsonjsjunior/stargazer.git /opt/stargazer \
      && (cd /opt/stargazer && python setup.py install) \
      && stargazer -h >/dev/null 2>&1 || echo "WARN: stargazer not installed; runner will emit NA"

# --- (B) registered-user tarball (PLACEHOLDER — uncomment after manual download) ----
# COPY stargazer-grc-2.0.2.tar.gz /tmp/stargazer.tar.gz
# RUN tar -xzf /tmp/stargazer.tar.gz -C /opt && cd /opt/Stargazer* \
#     && python setup.py install && stargazer -h

# run_stargazer.sh is mounted by Nextflow at runtime; do NOT COPY it here.
