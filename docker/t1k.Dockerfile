# T1K benchmark image  ->  pgxbench/t1k:1.0.5
#   docker build -f benchmark/docker/t1k.Dockerfile -t pgxbench/t1k:1.0.5 .
#
# Repo/docs: https://github.com/mourisl/T1K  (tag v1.0.5)
#            https://github.com/mourisl/T1K/blob/v1.0.5/README.md
# T1K genotypes highly-polymorphic genes; for CYP2D6 it is run against a PharmVar
# allele-definition reference (alleles named "CYP2D6*ABC").
#   Build ref: perl t1k-build.pl -f cyp2d6_pharmvar.fa -o /opt/t1kref --prefix cyp2d6
#              -> /opt/t1kref/cyp2d6_dna_seq.fa  (+ _rna_seq.fa, coord files if -g given)
#   Genotype : run-t1k (-1/-2 fastq | -b bam -c coord) -f <seq.fa> -o P --od D --preset kir-wgs
#   Output   : <P>_genotype.tsv  cols: gene num_alleles allele1 abund1 qual1 allele2 ...
#              (quality<=0 alleles should be ignored).
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential zlib1g-dev libbz2-dev liblzma-dev \
        git perl ca-certificates curl unzip \
        samtools bcftools tabix \
    && rm -rf /var/lib/apt/lists/*

# Compile T1K from the pinned v1.0.5 tag (depends on pthreads + zlib).
RUN git clone --depth 1 --branch v1.0.5 https://github.com/mourisl/T1K.git /opt/T1K \
    && make -C /opt/T1K \
    && ln -s /opt/T1K/run-t1k /usr/local/bin/run-t1k \
    && ln -s /opt/T1K/t1k /usr/local/bin/t1k \
    && ln -s /opt/T1K/t1k-build.pl /usr/local/bin/t1k-build.pl
ENV PATH=/opt/T1K:$PATH

# --- Build the CYP2D6 PharmVar reference into the image -----------------------------
# PharmVar gene ZIPs require an API key, so we cannot auto-download the latest set here.
# Two supported paths:
#   (A) [default] derive the allele-definition FASTA from the PharmVar GRCh38 reference
#       sequence + the per-allele variant VCFs that the build context provides, OR
#   (B) drop a prebuilt PharmVar CYP2D6 allele FASTA (alleles named "CYP2D6*ABC") into
#       the build context as cyp2d6_pharmvar.fa and uncomment the COPY below.
# If the reference cannot be built, the runner emits NA (it checks for the seq file).
#
# COPY cyp2d6_pharmvar.fa /opt/t1ksrc/cyp2d6_pharmvar.fa
RUN mkdir -p /opt/t1kref \
    && if [ -s /opt/t1ksrc/cyp2d6_pharmvar.fa ]; then \
         perl /opt/T1K/t1k-build.pl -f /opt/t1ksrc/cyp2d6_pharmvar.fa \
              -o /opt/t1kref --prefix cyp2d6 ; \
       else \
         echo "WARN: no PharmVar CYP2D6 FASTA in build context; T1K runner will emit NA" ; \
       fi

ENV T1K_REF_DIR=/opt/t1kref T1K_GENE=cyp2d6
# run_t1k.sh is mounted by Nextflow at runtime; do NOT COPY it here.
