# Cyrius benchmark image  ->  pgxbench/cyrius:1.1.1
#   docker build -f benchmark/docker/cyrius.Dockerfile -t pgxbench/cyrius:1.1.1 .
#
# Docs/repo: https://github.com/Illumina/Cyrius  (tag v1.1.1)
#            https://github.com/Illumina/Cyrius/blob/master/README.md
# Cyrius genotypes CYP2D6 ONLY, from a whole-genome BAM/CRAM:
#   python3 star_caller.py --manifest M --genome [19|37|38] \
#           --prefix P --outDir D [--reference ref.fa] [--threads N]
#   * manifest = text file, one absolute BAM/CRAM path per line.
#   * output  = <outDir>/<prefix>.tsv with columns  Sample <TAB> Genotype <TAB> Filter
#     Genotype e.g. "*1/*4", "*3/*68+*4", "*1x2/*2"; "None" = no-call.
FROM python:3.8-slim

# build toolchain for C-extension pip wheels (pytabix/mappy/pysam)
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential gcc g++ make zlib1g-dev libbz2-dev liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y --no-install-recommends \
        samtools tabix git \
    && rm -rf /var/lib/apt/lists/*

# Cyrius is distributed as a repo (star_caller.py + data/). Pin the v1.1.1 tag.
# numpy/scipy/pysam/statsmodels are the runtime deps used by the caller.
RUN pip install --no-cache-dir "numpy<2" scipy pysam statsmodels \
    && git clone --depth 1 --branch v1.1.1 https://github.com/Illumina/Cyrius.git /opt/Cyrius \
    && python3 /opt/Cyrius/star_caller.py --help >/dev/null

ENV CYRIUS_DIR=/opt/Cyrius
# run_cyrius.sh is mounted by Nextflow at runtime; do NOT COPY it here.
