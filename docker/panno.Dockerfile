# PAnno benchmark image  ->  pgxbench/panno:0.2.0
#   docker build -f benchmark/docker/panno.Dockerfile -t pgxbench/panno:0.2.0 .
#
# Docs/repo: https://github.com/PreMedKB/PAnno
#            https://pypi.org/project/panno/
#            https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2023.1008330/full
# PAnno parses a germline GRCh38 VCF + a biogeographic population and infers diplotypes
# (52 genes) + CPIC/DPWG phenotypes, rendering an HTML report:
#   panno -s <sample> -i <germline_vcf> -p <POP> -o <outdir>   ->  <outdir>/<sample>.html
#   population is one of: AAC AME EAS EUR LAT NEA OCE SAS SSA  (we default to EUR).
# PAnno writes ONLY the HTML report (diplotypes are computed in-memory then rendered),
# so the runner parses the per-gene diplotype out of the rendered HTML table.
# GRCh38 only (built-in definitions are GRCh38).
#
# The image tag is 0.2.0 per the harness convention, but we pin the PyPI package to a
# concrete release (0.3.1 is the current published version with the documented CLI).
FROM python:3.9-slim

# build toolchain for C-extension pip wheels (pytabix/mappy/pysam)
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential gcc g++ make zlib1g-dev libbz2-dev liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y --no-install-recommends \
        tabix bcftools \
    && rm -rf /var/lib/apt/lists/*

# pip install panno per the README (PyPI: panno). Pin to a concrete release.
RUN pip install --no-cache-dir "panno==0.3.1" \
 && python3 -c "import panno"            # smoke-check the package imports

# run_panno.sh is mounted by Nextflow at runtime; do NOT COPY it here.
