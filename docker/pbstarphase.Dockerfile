# pb-StarPhase benchmark image.
#   https://github.com/PacificBiosciences/pb-StarPhase
#   https://github.com/PacificBiosciences/pb-StarPhase/blob/main/docs/user_guide.md
#   bioconda: https://bioconda.github.io/recipes/pbstarphase/README.html
#
# A phase-aware pharmacogenomic diplotyper for PacBio HiFi data: ~21 CPIC genes
# (incl. CYP2D6 SV, HLA-A/B). Needs a pre-built allele database (shipped in the
# repo data/ folder) + a reference FASTA + a (bgzipped, .tbi-indexed) VCF; the BAM
# is required for the SV/HLA/CYP2D6 callers. The runner installs samtools/bcftools
# for slicing and produces the VCF the diplotyper expects.
#
#   docker build -f benchmark/docker/pbstarphase.Dockerfile -t pgxbench/pbstarphase:1.0.0 .
#
# Tool version is pinned via bioconda; the bundled DB pins the allele definitions.
FROM mambaorg/micromamba:1.5.8

# Tool itself is pbstarphase 2.x on bioconda; the image tag (1.0.0) is the
# benchmark-harness tag, independent of the upstream tool version.
ARG PBSTARPHASE_VERSION=2.1.0
# Pre-built database shipped in the repo (PharmVar/CPIC/HLA snapshot).
ARG PBSTARPHASE_DB_URL=https://raw.githubusercontent.com/PacificBiosciences/pb-StarPhase/main/data/v2.1.0/pbstarphase_20260622.json.gz

USER root
RUN micromamba install -y -n base -c bioconda -c conda-forge \
        pbstarphase=${PBSTARPHASE_VERSION} \
        samtools bcftools tabix curl \
    && micromamba clean -a -y

ENV PATH=/opt/conda/bin:$PATH

# Ship the pre-built allele database inside the image so the runner is offline.
ENV PBSTARPHASE_DB=/opt/pbstarphase/db.json.gz
RUN mkdir -p /opt/pbstarphase \
    && curl -fsSL "${PBSTARPHASE_DB_URL}" -o "${PBSTARPHASE_DB}" \
    && pbstarphase --version
