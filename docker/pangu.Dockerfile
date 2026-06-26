# pangu benchmark image.
#   https://github.com/PacificBiosciences/pangu
#   bioconda: https://bioconda.github.io/recipes/pangu/README.html
#
# PacBio star-typer for CYP2D6 ONLY (incl. structural/CNV alleles) from a HiFi BAM
# aligned to GRCh38. Emits a <prefix>_report.json whose [0].diplotype is e.g.
# "CYP2D6 *4x2/*5". No reference FASTA required for BAM input; samtools is bundled
# so the runner can index/validate the slice.
#
#   docker build -f benchmark/docker/pangu.Dockerfile -t pgxbench/pangu:0.2.1 .
FROM mambaorg/micromamba:1.5.8

# pangu upstream version pinned to 0.2.x; image tag matches the harness pin.
ARG PANGU_VERSION=0.2.8

USER root
RUN micromamba install -y -n base -c bioconda -c conda-forge \
        pangu=${PANGU_VERSION} \
        samtools \
    && micromamba clean -a -y

ENV PATH=/opt/conda/bin:$PATH
RUN pangu --help >/dev/null 2>&1 || pangu -h >/dev/null 2>&1 || true
