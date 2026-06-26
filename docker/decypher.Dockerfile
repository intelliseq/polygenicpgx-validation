# deCYPher benchmark image.
#   Preprint: https://www.biorxiv.org/content/10.1101/2025.10.13.681303v2.full
#   DOI: 10.1101/2025.10.13.681303
#
# STATUS (2026-06): deCYPher is a 2025 bioRxiv preprint. A star-allele-resolution
# framework that operates on HAPLOTYPE-RESOLVED long-read ASSEMBLIES (phased
# contigs / FASTA), not on a per-gene aligned BAM, for the PharmVar 1A genes
# (CYP2B6/2C9/2C19/2D6/3A5/4F2, DPYD, NUDT15, SLCO1B1). As of this build no public
# GitHub repository or released CLI/binary was discoverable (searched github +
# bioRxiv code-availability). This image therefore installs the *documented*
# toolchain the method depends on (minimap2/samtools/python) so the runner can be
# wired up the moment a CLI is published; until then run_decypher.sh emits NA
# because (a) the decypher binary is absent and (b) the bundle carries an aligned
# BAM slice, not the phased-assembly FASTA deCYPher requires.
#
#   docker build -f benchmark/docker/decypher.Dockerfile -t pgxbench/decypher:0.1.0 .
#
# TODO(verify): replace the install block below with the upstream install once the
# repo is released (expected something like `pip install decypher` or a clone +
# `pip install .`), and point DECYPHER_BIN at the real entrypoint.
FROM mambaorg/micromamba:1.5.8

USER root
RUN micromamba install -y -n base -c bioconda -c conda-forge \
        python=3.10 minimap2 samtools bcftools \
    && micromamba clean -a -y

ENV PATH=/opt/conda/bin:$PATH

# Entrypoint name the runner probes for; unset/absent today -> runner emits NA.
ENV DECYPHER_BIN=decypher
