# Chinook / EPI2ME wf-pgx benchmark image.
#   Workflow index: https://epi2me.nanoporetech.com/wfindex/
#   ONT PGx protocol: https://nanoporetech.com/document/pgx-sequencing-workflow-with-adaptive-sampling-blood-cells-saliva-
#   CYP2D6 resolution write-up: https://nanoporetech.com/resource-centre/full-cyp2d6-and-pharmacogene-resolution-with-oxford-nanopore-long-reads
#
# STATUS (2026-06): wf-pgx is an Oxford Nanopore EPI2ME Nextflow workflow. It calls
# SNV/indels and runs PharmCAT for most PGx genes, and uses "Chinook" — a targeted
# CYP2D6 assembly/SV caller — for CYP2D6. As of this build there is NO public
# `epi2me-labs/wf-pgx` GitHub repo (github.com/epi2me-labs/wf-pgx -> 404) and the
# Chinook binary is not separately distributed; the workflow is delivered through
# the EPI2ME desktop/Nextflow channel. Running nextflow-in-nextflow is also out of
# scope for this harness. This image therefore ships nextflow + the ONT toolchain so
# the workflow can be slotted in once published, but run_chinook.sh emits NA today
# (no wf-pgx/Chinook entrypoint present).
#
#   docker build -f benchmark/docker/chinook.Dockerfile -t pgxbench/wf-pgx:latest .
#
# TODO(verify): once epi2me-labs/wf-pgx is public, either (a) `nextflow run
# epi2me-labs/wf-pgx --bam ... --ref ...` against the gene BAM, or (b) invoke the
# standalone Chinook CYP2D6 caller; point CHINOOK_BIN / WFPGX_DIR at it.
FROM mambaorg/micromamba:1.5.8

USER root
RUN micromamba install -y -n base -c bioconda -c conda-forge \
        nextflow samtools minimap2 \
    && micromamba clean -a -y

ENV PATH=/opt/conda/bin:$PATH

# Entrypoints the runner probes for; absent today -> runner emits NA.
ENV CHINOOK_BIN=chinook
ENV WFPGX_DIR=/opt/wf-pgx
