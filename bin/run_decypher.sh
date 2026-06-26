#!/usr/bin/env bash
# Runner contract (shared by every run_<tool>.sh):
#   args: <sample> <gene> <build> <bundle_dir> <out_raw>
#   writes <out_raw> as a 2-col TSV with keys: tool sample gene build status diplotype phenotype
#   status: OK (a call was produced) | NA (tool/build/gene unsupported or input missing) | ERROR
#
# deCYPher: bioRxiv 2025 (DOI 10.1101/2025.10.13.681303)
#   https://www.biorxiv.org/content/10.1101/2025.10.13.681303v2.full
#   Star-allele resolution for the PharmVar 1A genes (CYP2B6/2C9/2C19/2D6/3A5/4F2,
#   DPYD, NUDT15, SLCO1B1) from HAPLOTYPE-RESOLVED long-read ASSEMBLIES (phased
#   contigs), NOT from an aligned per-gene BAM.
#
# STATUS (2026-06): no public GitHub repo / released CLI was discoverable. This
# runner is written faithfully against what the preprint documents:
#   - it needs a phased-assembly FASTA (a `decypher` binary consuming contigs), and
#   - the benchmark bundle only carries an aligned BAM slice, not phased contigs.
# Because the binary is absent AND the assembly input is absent, the runner emits NA
# (never invents a call). Replace the probe + invocation below once a CLI is
# released. The honest design keeps the harness contract satisfied today.
set -uo pipefail
sample=$1 gene=$2 build=$3 bundle=$4 out=$5
TOOL=decypher

emit() { # status diplotype [phenotype]
  printf 'tool\t%s\nsample\t%s\ngene\t%s\nbuild\t%s\nstatus\t%s\ndiplotype\t%s\nphenotype\t%s\n' \
    "$TOOL" "$sample" "$gene" "$build" "$1" "${2:-}" "${3:-}" > "$out"
}

# PharmVar 1A genes deCYPher targets (cyp4f2 not in our panel; vkorc1/tpmt unsupported).
case "$gene" in
  cyp2d6|cyp2c19|cyp2c9|cyp3a5|cyp2b6|slco1b1|nudt15|dpyd) ;;
  *) emit NA ""; exit 0 ;;
esac
# GRCh38-only.
[ "$build" = grch38 ] || { emit NA ""; exit 0; }

bam="$bundle/${sample}.${gene}.longread.bam"
# Long-read slice only exists for the 3-sample sub-cohort; absent -> NA.
[ -s "$bam" ] || { emit NA ""; exit 0; }

# deCYPher CLI absent (no released binary) -> honest NA.
DECYPHER_BIN="${DECYPHER_BIN:-decypher}"
command -v "$DECYPHER_BIN" >/dev/null 2>&1 || { emit NA ""; exit 0; }

# deCYPher requires a phased-assembly FASTA (haplotype-resolved contigs). The bundle
# does not carry one, so even with a binary present we have no valid input -> NA.
asm="$bundle/${sample}.${gene}.asm.fasta"
[ -s "$asm" ] || { emit NA ""; exit 0; }

# --- Below is the intended invocation once a CLI/assembly input exist (untested). ---
# work=$(mktemp -d); ref="${PGXBENCH_GRCH38_REF:-}"
# "$DECYPHER_BIN" --gene "$gene" --assembly "$asm" --reference "$ref" \
#     --out "$work/out.json" >/dev/null 2>&1 || { emit ERROR ""; rm -rf "$work"; exit 0; }
# diplo=$(python3 -c '...parse out.json...')
# rm -rf "$work"
# [ -n "$diplo" ] && emit OK "$diplo" "" || emit ERROR ""

emit NA ""
exit 0
