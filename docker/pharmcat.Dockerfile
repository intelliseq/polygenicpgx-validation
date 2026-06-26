# PharmCAT benchmark image  ->  pgxbench/pharmcat:2.15.5
#   docker build -f benchmark/docker/pharmcat.Dockerfile -t pgxbench/pharmcat:2.15.5 .
#
# Docs: https://pharmcat.clinpgx.org/using/Running-PharmCAT-Pipeline/
# The official PharmGKB image already ships the JRE, the PharmCAT jar, and the
# `pharmcat_pipeline` wrapper (preprocessor + named-allele matcher + phenotyper +
# reporter) on PATH, plus bcftools/bgzip/tabix that the preprocessor needs.
FROM pgkb/pharmcat:2.15.5

# samtools is handy for any VCF housekeeping in the runner; the base image already
# carries bcftools + htslib (bgzip/tabix) used by the preprocessor.
USER root
RUN (apt-get update && apt-get install -y --no-install-recommends samtools \
     && rm -rf /var/lib/apt/lists/*) || true

# run_pharmcat.sh is mounted by Nextflow at runtime; do NOT COPY it here.
