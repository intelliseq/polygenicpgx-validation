# StellarPGx benchmark image  ->  pgxbench/stellarpgx:1.2.7
#   docker build -f benchmark/docker/stellarpgx.Dockerfile -t pgxbench/stellarpgx:1.2.7 .
#
# Repo/docs: https://github.com/SBIMB/StellarPGx
#            https://github.com/SBIMB/StellarPGx/blob/master/README.md
# StellarPGx is a Nextflow pipeline that calls CYP450 star alleles from
# *high-coverage whole-genome* BAM/CRAM (CYP2D6/2C19/2C9/3A5/2B6 + others).
#   nextflow run main.nf -profile standard --build [hg38|hg19|b37] --gene cyp2d6 \
#           --in_bam "/data/<sample>*{bam,bai}" --ref_file REF --out_dir D --format compressed
#   Result: <out_dir>/<gene>/alleles/<sample>_<gene>.alleles  (text report with a
#   "Result:" diplotype line, e.g. "*1/*4", plus "Metaboliser status:" line).
#
# IMPORTANT: StellarPGx requires WHOLE-GENOME input (it derives CNV from read-depth
# ratios against control regions); it cannot run on a gene-slice CRAM. The runner is
# honest about this: it only invokes the pipeline when a whole-genome BAM is supplied,
# otherwise it emits status=NA. The image bundles the runtime (Nextflow + the pinned
# StellarPGx repo + caller deps) so it CAN run if whole-genome input is ever provided.
FROM nextflow/nextflow:23.10.1

# Tooling the pipeline's processes shell out to (graphtyper is vendored by the
# SBIMB container image; here we provide the standard read/variant utilities and a
# Python 3.8 the custom caller modules need).
RUN if command -v microdnf >/dev/null 2>&1; then \
        microdnf install -y git python38 procps-ng tar gzip which findutils \
          && microdnf clean all ; \
    else \
        (yum install -y git python3 procps-ng tar gzip which findutils || true) ; \
    fi

# Pin the StellarPGx repo. Upstream nextflow.config manifest reads 1.2.8 in master;
# we pin the 1.2.7 release tag to match the pgxbench/stellarpgx:1.2.7 image tag.
# (The SBIMB/Docker-Hub image twesigomwedavid/stellarpgx-dev provides graphtyper +
# the genome-graph deps; the runner can point the pipeline's container at it.)
RUN git clone --depth 1 --branch 1.2.7 https://github.com/SBIMB/StellarPGx.git /opt/StellarPGx \
      || git clone --depth 1 https://github.com/SBIMB/StellarPGx.git /opt/StellarPGx

ENV STELLARPGX_DIR=/opt/StellarPGx \
    STELLARPGX_IMAGE=twesigomwedavid/stellarpgx-dev:latest
# run_stellarpgx.sh is mounted by Nextflow at runtime; do NOT COPY it here.
