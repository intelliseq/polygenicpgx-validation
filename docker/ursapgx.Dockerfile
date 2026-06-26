# ursaPGx benchmark image  ->  pgxbench/ursapgx:1.0.0
#   docker build -f benchmark/docker/ursapgx.Dockerfile -t pgxbench/ursapgx:1.0.0 .
#
# Docs/repo: https://github.com/coriell-research/ursaPGx
#            https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2024.1351620/full
# ursaPGx is an R/Bioconductor package that assigns PHASED diplotype calls from a
# phased single- or multi-sample indexed VCF using PharmVar star-allele definitions.
# It is GRCh38/GRCh37 capable, but our bundle only carries PHASED GRCh38 slices, so the
# runner restricts it to grch38. R API (callDiplotypes wraps the pipeline):
#   library(ursaPGx)
#   res <- callDiplotypes("phased.vcf.gz", gene = "CYP2C19", phased = TRUE)
#   # res is a DataFrame: rows = samples, one column = the gene, value e.g. "*1|*2"
# CYP2D6 in ursaPGx is handled via Cyrius (BAM/CRAM) not VCF, so the driver reports it
# uncallable from a VCF slice -> NA.
#
# Pin a Bioconductor release for reproducibility (RELEASE_3_18 ~ R 4.3 / Bioc 3.18).
FROM bioconductor/bioconductor_docker:RELEASE_3_18

RUN apt-get update && apt-get install -y --no-install-recommends \
        tabix bcftools \
    && rm -rf /var/lib/apt/lists/*

# VariantAnnotation is the Bioconductor backbone ursaPGx inherits from; install it from
# Bioconductor, then install ursaPGx from GitHub (pinned commit-less main via remotes).
RUN R -q -e 'install.packages("remotes", repos="https://cloud.r-project.org")' \
 && R -q -e 'BiocManager::install(c("VariantAnnotation","GenomicRanges","Biostrings","S4Vectors"), ask=FALSE, update=FALSE)' \
 && R -q -e 'remotes::install_github("coriell-research/ursaPGx", upgrade="never")' \
 && R -q -e 'library(ursaPGx)'   # smoke-check the package loads

# Tool-specific R driver (NOT a bin/ runner) baked into the image.
COPY ursapgx_call.R /opt/ursapgx_call.R

ENV URSAPGX_DRIVER=/opt/ursapgx_call.R
# run_ursapgx.sh is mounted by Nextflow at runtime; do NOT COPY it here.
