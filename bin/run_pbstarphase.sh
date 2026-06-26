#!/usr/bin/env bash
# Runner contract (shared by every run_<tool>.sh):
#   args: <sample> <gene> <build> <bundle_dir> <out_raw>
#   writes <out_raw> as a 2-col TSV with keys: tool sample gene build status diplotype phenotype
#   status: OK (a call was produced) | NA (tool/build/gene unsupported or input missing) | ERROR
#
# pb-StarPhase: https://github.com/PacificBiosciences/pb-StarPhase
#   https://github.com/PacificBiosciences/pb-StarPhase/blob/main/docs/user_guide.md
#   Phase-aware diplotyper for PacBio HiFi. ~21 CPIC genes incl. CYP2D6 SV.
#     pbstarphase diplotype --database DB --reference REF --vcf VCF[.tbi] \
#                  [--bam BAM] --output-calls OUT.json
#   Output JSON: gene_details.<GENE>.diplotypes[0].diplotype  == "*1/*4" etc.
#   It only assigns diplotypes (no phenotype) so we leave phenotype blank.
# GRCh38-only; needs the per-gene long-read BAM slice (3-sample sub-cohort only).
set -uo pipefail
sample=$1 gene=$2 build=$3 bundle=$4 out=$5
TOOL=pbstarphase

emit() { # status diplotype [phenotype]
  printf 'tool\t%s\nsample\t%s\ngene\t%s\nbuild\t%s\nstatus\t%s\ndiplotype\t%s\nphenotype\t%s\n' \
    "$TOOL" "$sample" "$gene" "$build" "$1" "${2:-}" "${3:-}" > "$out"
}

# Map our lowercase gene id -> pb-StarPhase (uppercase) gene_details key.
declare -A GMAP=(
  [cyp2d6]=CYP2D6 [cyp2c19]=CYP2C19 [cyp2c9]=CYP2C9 [cyp3a5]=CYP3A5
  [cyp2b6]=CYP2B6 [slco1b1]=SLCO1B1 [vkorc1]=VKORC1 [tpmt]=TPMT
  [nudt15]=NUDT15 [dpyd]=DPYD
)
G=${GMAP[$gene]:-}
[ -n "$G" ] || { emit NA ""; exit 0; }

# GRCh38-only tool.
[ "$build" = grch38 ] || { emit NA ""; exit 0; }

bam="$bundle/${sample}.${gene}.longread.bam"
ref="${PGXBENCH_GRCH38_REF:-}"
db="${PBSTARPHASE_DB:-/opt/pbstarphase/db.json.gz}"

# Long-read slice only exists for the 3-sample sub-cohort; absent -> NA.
[ -s "$bam" ] || { emit NA ""; exit 0; }
[ -n "$ref" ] && [ -s "$ref" ] || { emit NA ""; exit 0; }
[ -s "$db" ] || { emit ERROR ""; exit 0; }

work=$(mktemp -d)
calls="$work/calls.json"

# pb-StarPhase wants a (bgzipped, .tbi-indexed) VCF for variant-based genes plus the
# BAM for the SV/CYP2D6 callers. Derive a small VCF from the gene BAM slice.
[ -s "${bam}.bai" ] || samtools index "$bam" >/dev/null 2>&1
vcf="$work/calls.vcf.gz"
if bcftools mpileup -f "$ref" "$bam" 2>/dev/null \
     | bcftools call -mv -Oz -o "$vcf" 2>/dev/null \
   && tabix -p vcf "$vcf" >/dev/null 2>&1; then
  VCFARG=(--vcf "$vcf")
else
  VCFARG=()   # variant VCF optional in recent pb-StarPhase; BAM still drives SV/CYP2D6
fi

pbstarphase diplotype \
  --database "$db" --reference "$ref" \
  "${VCFARG[@]}" --bam "$bam" \
  --output-calls "$calls" >/dev/null 2>&1 \
  || { emit ERROR ""; rm -rf "$work"; exit 0; }

[ -s "$calls" ] || { emit ERROR ""; rm -rf "$work"; exit 0; }

diplo=$(python3 - "$calls" "$G" <<'PY'
import sys, json
diplo = ""
try:
    d = json.load(open(sys.argv[1]))
    g = d.get("gene_details", {}).get(sys.argv[2], {})
    dips = g.get("diplotypes") or g.get("simple_diplotypes") or []
    if dips:
        first = dips[0]
        diplo = first.get("diplotype") or ""
        if not diplo and first.get("hap1") is not None:
            diplo = f"{first.get('hap1','')}/{first.get('hap2','')}"
except Exception:
    pass
print((diplo or "").replace(" ", "_"))
PY
)
rm -rf "$work"

[ -n "$diplo" ] && emit OK "$diplo" "" || emit ERROR ""
