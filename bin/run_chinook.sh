#!/usr/bin/env bash
# Runner contract (shared by every run_<tool>.sh):
#   args: <sample> <gene> <build> <bundle_dir> <out_raw>
#   writes <out_raw> as a 2-col TSV with keys: tool sample gene build status diplotype phenotype
#   status: OK (a call was produced) | NA (tool/build/gene unsupported or input missing) | ERROR
#
# Chinook / EPI2ME wf-pgx: Oxford Nanopore PGx workflow.
#   https://epi2me.nanoporetech.com/wfindex/
#   https://nanoporetech.com/resource-centre/full-cyp2d6-and-pharmacogene-resolution-with-oxford-nanopore-long-reads
#   wf-pgx calls SNV/indels + PharmCAT for most genes; "Chinook" is its targeted
#   CYP2D6 assembly/SV caller. We focus on CYP2D6 (other genes -> NA).
#
# STATUS (2026-06): github.com/epi2me-labs/wf-pgx is NOT public (404) and Chinook is
# not separately distributed; the workflow ships via the EPI2ME channel only.
# Running a Nextflow workflow inside this harness is also out of scope. So this
# runner probes for a Chinook binary / wf-pgx checkout and, finding neither, emits
# NA (never invents a call). The intended CYP2D6 invocation is sketched below for
# when the workflow becomes runnable.
set -uo pipefail
sample=$1 gene=$2 build=$3 bundle=$4 out=$5
TOOL=chinook

emit() { # status diplotype [phenotype]
  printf 'tool\t%s\nsample\t%s\ngene\t%s\nbuild\t%s\nstatus\t%s\ndiplotype\t%s\nphenotype\t%s\n' \
    "$TOOL" "$sample" "$gene" "$build" "$1" "${2:-}" "${3:-}" > "$out"
}

# Chinook is the CYP2D6 caller; this wrapper is scoped to CYP2D6.
[ "$gene" = cyp2d6 ] || { emit NA ""; exit 0; }
# GRCh38-only (ONT long reads aligned to GRCh38).
[ "$build" = grch38 ] || { emit NA ""; exit 0; }

bam="$bundle/${sample}.${gene}.longread.bam"
# Long-read slice only exists for the 3-sample sub-cohort; absent -> NA.
[ -s "$bam" ] || { emit NA ""; exit 0; }

# Chinook binary / wf-pgx checkout absent today -> honest NA.
CHINOOK_BIN="${CHINOOK_BIN:-chinook}"
if ! command -v "$CHINOOK_BIN" >/dev/null 2>&1 \
   && [ ! -d "${WFPGX_DIR:-/opt/wf-pgx}/.git" ] \
   && [ ! -f "${WFPGX_DIR:-/opt/wf-pgx}/main.nf" ]; then
  emit NA ""; exit 0
fi

# --- Intended invocation once Chinook/wf-pgx are runnable (untested). ---
# work=$(mktemp -d); ref="${PGXBENCH_GRCH38_REF:-}"
# [ -s "${bam}.bai" ] || samtools index "$bam" >/dev/null 2>&1
# "$CHINOOK_BIN" --bam "$bam" --reference "$ref" --out "$work/cyp2d6.json" \
#     >/dev/null 2>&1 || { emit ERROR ""; rm -rf "$work"; exit 0; }
# diplo=$(python3 -c '...parse Chinook CYP2D6 diplotype...')
# rm -rf "$work"
# [ -n "$diplo" ] && emit OK "$diplo" "" || emit ERROR ""

emit NA ""
exit 0
