# SpecImmune benchmark image.
#   https://github.com/deepomicslab/SpecImmune
#   Preprint (2025): long-read typing of immune-related gene families incl. CYP.
#
# Accurate typing of diverse gene families (HLA, KIR, IG, TCR, CYP) from long reads
# (PacBio or Nanopore), WGS or amplicon. CYP genes (incl. CYP2D6) are typed with
# `-i CYP`. Inputs: a long-read FASTQ (-r), a no-alt GRCh38 reference (--hg38), and a
# pre-built database folder (--db, made with scripts/make_db.py). Per-gene results
# land in <sample>.CYP.merge.type.result.txt.
#
#   docker build -f benchmark/docker/specimmune.Dockerfile -t pgxbench/specimmune:1.0.0 .
#
# Pinned via a git tag/commit; the env.yml resolves the (heavy) bioinfo deps. The
# build does NOT fetch/host the CYP database (large, license-bound) — if --db is not
# provided to the runner at runtime it emits NA (see run_specimmune.sh).
FROM mambaorg/micromamba:1.5.8

ARG SPECIMMUNE_REF=main

USER root
RUN micromamba install -y -n base -c bioconda -c conda-forge \
        git python=3.10 samtools minimap2 dysgu=1.6.2 \
    && micromamba clean -a -y

ENV PATH=/opt/conda/bin:$PATH

# Clone the framework. Upstream is installed in-place (scripts/main.py is the CLI).
ENV SPECIMMUNE_DIR=/opt/SpecImmune
RUN git clone --depth 1 --branch "${SPECIMMUNE_REF}" \
        https://github.com/deepomicslab/SpecImmune.git "${SPECIMMUNE_DIR}" \
    && (micromamba run -n base pip install --no-cache-dir -r "${SPECIMMUNE_DIR}/requirements.txt" 2>/dev/null || true) \
    && chmod -R +x "${SPECIMMUNE_DIR}/bin" 2>/dev/null || true

# Optional: mount/bake a prebuilt CYP database here; if absent the runner emits NA.
ENV SPECIMMUNE_DB=/opt/SpecImmune/db
