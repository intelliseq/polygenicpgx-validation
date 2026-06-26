#!/usr/bin/env bash
# Runner contract (shared by every run_<tool>.sh):
#   args: <sample> <gene> <build> <bundle_dir> <out_raw>
#   writes <out_raw> as a 2-col TSV with keys: tool sample gene build status diplotype phenotype
#   status: OK (a call was produced) | NA (tool/build/gene unsupported or input missing) | ERROR
#
# SpecImmune: https://github.com/deepomicslab/SpecImmune  (bioRxiv 2025)
#   Long-read typing of immune + CYP gene families (PacBio/Nanopore, WGS/amplicon).
#   CYP genes (incl. CYP2D6) are typed with `-i CYP`:
#     python3 scripts/main.py -r <fastq> -j <threads> -i CYP -n <sample> \
#         -o <outdir> --hg38 <no_alt_ref> --db <db_dir> -y <pacbio|nanopore>
#   Per-gene results: <outdir>/<sample>.CYP.merge.type.result.txt
#   SpecImmune types per CYP locus; we extract the requested gene's allele row.
# GRCh38-only; needs the long-read slice (3-sample sub-cohort) + a prebuilt CYP DB.
set -uo pipefail
sample=$1 gene=$2 build=$3 bundle=$4 out=$5
TOOL=specimmune

emit() { # status diplotype [phenotype]
  printf 'tool\t%s\nsample\t%s\ngene\t%s\nbuild\t%s\nstatus\t%s\ndiplotype\t%s\nphenotype\t%s\n' \
    "$TOOL" "$sample" "$gene" "$build" "$1" "${2:-}" "${3:-}" > "$out"
}

# SpecImmune CYP family -> our panel intersection (uppercase locus label it reports).
declare -A GMAP=(
  [cyp2d6]=CYP2D6 [cyp2c19]=CYP2C19 [cyp2c9]=CYP2C9 [cyp3a5]=CYP3A5
  [cyp2b6]=CYP2B6
  # SLCO1B1/VKORC1/TPMT/NUDT15/DPYD are not in the SpecImmune CYP gene family.
)
G=${GMAP[$gene]:-}
[ -n "$G" ] || { emit NA ""; exit 0; }

# GRCh38-only.
[ "$build" = grch38 ] || { emit NA ""; exit 0; }

bam="$bundle/${sample}.${gene}.longread.bam"
ref="${PGXBENCH_GRCH38_REF:-}"
db="${SPECIMMUNE_DB:-/opt/SpecImmune/db}"
SPECIMMUNE_DIR="${SPECIMMUNE_DIR:-/opt/SpecImmune}"
plat="${SPECIMMUNE_PLATFORM:-pacbio}"   # bundle long reads are PacBio HiFi by default

# Long-read slice only exists for the 3-sample sub-cohort; absent -> NA.
[ -s "$bam" ] || { emit NA ""; exit 0; }
[ -n "$ref" ] && [ -s "$ref" ] || { emit NA ""; exit 0; }
# Prebuilt CYP DB and the framework checkout are required; absent -> NA (not bundled).
[ -d "$db" ] || { emit NA ""; exit 0; }
[ -f "$SPECIMMUNE_DIR/scripts/main.py" ] || { emit NA ""; exit 0; }

work=$(mktemp -d)
# SpecImmune consumes a FASTQ (-r); derive reads from the BAM slice.
[ -s "${bam}.bai" ] || samtools index "$bam" >/dev/null 2>&1
fq="$work/${sample}.fastq.gz"
samtools fastq "$bam" 2>/dev/null | gzip -c > "$fq" 2>/dev/null
[ -s "$fq" ] || { emit ERROR ""; rm -rf "$work"; exit 0; }

python3 "$SPECIMMUNE_DIR/scripts/main.py" \
  -r "$fq" -j 1 -i CYP -n "$sample" -o "$work" \
  --hg38 "$ref" --db "$db" -y "$plat" >/dev/null 2>&1 \
  || { emit ERROR ""; rm -rf "$work"; exit 0; }

res="$work/${sample}.CYP.merge.type.result.txt"
[ -s "$res" ] || res=$(ls "$work"/*CYP*type*result*.txt 2>/dev/null | head -1)
[ -n "$res" ] && [ -s "$res" ] || { emit ERROR ""; rm -rf "$work"; exit 0; }

# Extract the requested locus row; SpecImmune reports a per-gene allele/diplotype.
diplo=$(python3 - "$res" "$G" <<'PY'
import sys, re
res, gene = sys.argv[1], sys.argv[2]
diplo = ""
try:
    for line in open(res):
        s = line.rstrip("\n")
        if not s or s.startswith("#"):
            continue
        cols = s.split("\t")
        # row keyed by the locus name in any column
        if any(c.strip().upper() == gene for c in cols):
            # join the *-bearing allele tokens into a diplotype
            alleles = [c.strip() for c in cols if "*" in c]
            if alleles:
                diplo = "/".join(alleles[:2]) if len(alleles) > 1 else alleles[0]
            break
    # normalise gene-prefixed tokens "CYP2D6*4" -> "*4", keep xN/+
    diplo = re.sub(r"[A-Za-z0-9]+\*", "*", diplo)
except Exception:
    pass
print((diplo or "").replace(" ", "_"))
PY
)
rm -rf "$work"

[ -n "$diplo" ] && emit OK "$diplo" "" || emit ERROR ""
