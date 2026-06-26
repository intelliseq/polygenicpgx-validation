#!/usr/bin/env bash
# Launch the PGx benchmark harness with all data-source env vars set.
# Extra args are passed straight to nextflow (e.g. --samples NA12878 --tools polygenic -resume).
# Run from the repo root:  bash benchmark/run_benchmark.sh [nextflow args...]
set -uo pipefail
cd "$(dirname "$0")"           # -> benchmark/
REF="$PWD/work/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"

export PGXBENCH_GRCH38_REF="$REF"
export ALDY_REFERENCE="$REF"
export PGXBENCH_SEQ_INDEX="$PWD/work/sequence.index"
# CRAM decoding cache for samtools (ONT long-read CRAM etc.)
export REF_PATH="$REF"
# Long-read sub-cohort (2 of 3 verified public GRCh38 alignments; HG01190 has none).
export PGXBENCH_LONGREAD_NA12878="https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/PacBio_SequelII_CCS_11kb/HG001_GRCh38/HG001_GRCh38.haplotag.RTG.trio.bam"
export PGXBENCH_LONGREAD_HG00731="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/HG00731.hg38.cram"

echo "[run] REF=$REF"
echo "[run] nextflow run main.nf -profile docker $*"
nextflow run main.nf -profile docker "$@"
