#!/usr/bin/env bash
# Runner contract (shared by every run_<tool>.sh):
#   args: <sample> <gene> <build> <bundle_dir> <out_raw>
#   writes <out_raw> as a 2-col TSV with keys: tool sample gene build status diplotype phenotype
#   status: OK (a call was produced) | NA (tool/build/gene unsupported or input missing) | ERROR
#
# StellarPGx: https://github.com/SBIMB/StellarPGx  (Nextflow pipeline, v1.2.7)
#   nextflow run main.nf -profile standard --build [hg38|hg19|b37] --gene <gene> \
#           --in_bam "<dir>/<sample>*{bam,bai}" --ref_file REF --out_dir D [--format compressed]
#   Result file: <out_dir>/<gene>/alleles/<sample>_<gene>.alleles  (text report; the
#   diplotype follows a "Result:" line, the metaboliser status follows
#   "Metaboliser status:").
#
# HONEST LIMITATION: StellarPGx needs *whole-genome* high-coverage BAM/CRAM — it
# derives copy number from read-depth ratios vs control regions, so it CANNOT run on
# the gene-slice CRAM the bundle ships ($bundle/<sample>-<gene>.grch38.cram). We only
# invoke the pipeline if a genuine whole-genome BAM is provided via env
# PGXBENCH_WGS_BAM (or $bundle/<sample>.wgs.bam). Absent that -> status=NA (we do NOT
# fake a call off the slice).
set -uo pipefail
sample=$1 gene=$2 build=$3 bundle=$4 out=$5
TOOL=stellarpgx

emit() { # status diplotype [phenotype]
  printf 'tool\t%s\nsample\t%s\ngene\t%s\nbuild\t%s\nstatus\t%s\ndiplotype\t%s\nphenotype\t%s\n' \
    "$TOOL" "$sample" "$gene" "$build" "$1" "${2:-}" "${3:-}" > "$out"
}

# StellarPGx-supported genes from our panel (gene id passed lowercase to --gene).
case "$gene" in
  cyp2d6|cyp2c19|cyp2c9|cyp3a5|cyp2b6|slco1b1|nudt15|tpmt) : ;;
  *) emit NA ""; exit 0 ;;   # vkorc1/dpyd not modelled by StellarPGx
esac

# Map our build label -> StellarPGx --build value (chr-prefixed assemblies).
case "$build" in
  grch38) SBUILD=hg38 ;;
  grch37) SBUILD=hg19 ;;
  *)      emit NA ""; exit 0 ;;
esac

# Whole-genome BAM is mandatory; a gene-slice CRAM is insufficient for CNV-aware calls.
wgs="${PGXBENCH_WGS_BAM:-}"
[ -n "$wgs" ] || wgs="$bundle/${sample}.wgs.bam"
if [ ! -s "$wgs" ]; then
  # Be honest: no whole-genome input -> StellarPGx cannot run on the slice. NA.
  emit NA ""; exit 0
fi

ref="${PGXBENCH_GRCH38_REF:-}"
[ -n "$ref" ] && [ -s "$ref" ] || { emit NA ""; exit 0; }

work=$(mktemp -d)
indir="$work/data"; mkdir -p "$indir"
# StellarPGx globs "<dir>/<sample>*{bam,bai}". Stage links named for the sample.
ln -sf "$wgs" "$indir/${sample}.bam"
if [ -s "${wgs}.bai" ]; then ln -sf "${wgs}.bai" "$indir/${sample}.bam.bai"
elif [ -s "${wgs%.bam}.bai" ]; then ln -sf "${wgs%.bam}.bai" "$indir/${sample}.bam.bai"; fi
outdir="$work/results"

( cd "${STELLARPGX_DIR:-/opt/StellarPGx}" && \
  nextflow run main.nf -profile standard \
    --build "$SBUILD" --gene "$gene" \
    --in_bam "$indir/${sample}*{bam,bai}" \
    --ref_file "$ref" --out_dir "$outdir" \
    -with-docker "${STELLARPGX_IMAGE:-twesigomwedavid/stellarpgx-dev:latest}" \
) >/dev/null 2>&1 || { emit ERROR ""; rm -rf "$work"; exit 0; }

res=$(ls "$outdir/${gene}/alleles/"*".alleles" 2>/dev/null | head -1)
[ -n "$res" ] && [ -s "$res" ] || { emit ERROR ""; rm -rf "$work"; exit 0; }

read -r diplo pheno < <(python3 - "$res" <<'PY'
import sys, re
diplo = pheno = ""
lines = [l.rstrip("\n") for l in open(sys.argv[1])]
for i, l in enumerate(lines):
    s = l.strip()
    if s.lower().startswith("result:"):
        rest = s.split(":", 1)[1].strip()
        diplo = rest if rest else (lines[i+1].strip() if i+1 < len(lines) else "")
    elif s.lower().startswith("metaboliser status:") or s.lower().startswith("metabolizer status:"):
        rest = s.split(":", 1)[1].strip()
        pheno = rest if rest else (lines[i+1].strip() if i+1 < len(lines) else "")
# Keep only a star-allele diplotype (e.g. *1/*4, *5/*5, *1x2/*2); drop prose results.
if not re.match(r"^\*?\w+(/|\+|x).*$|^\*\d", diplo) and "*" not in diplo:
    diplo = ""
if diplo.lower() in ("indeterminate", "none", "na"):
    diplo = ""
print(diplo.replace(" ", "_"), pheno.replace(" ", "_"))
PY
)
rm -rf "$work"

[ -n "$diplo" ] && emit OK "$diplo" "$pheno" || emit ERROR ""
