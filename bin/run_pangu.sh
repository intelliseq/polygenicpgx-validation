#!/usr/bin/env bash
# Runner contract (shared by every run_<tool>.sh):
#   args: <sample> <gene> <build> <bundle_dir> <out_raw>
#   writes <out_raw> as a 2-col TSV with keys: tool sample gene build status diplotype phenotype
#   status: OK (a call was produced) | NA (tool/build/gene unsupported or input missing) | ERROR
#
# pangu: https://github.com/PacificBiosciences/pangu
#   CYP2D6 ONLY (other genes -> NA), from a PacBio HiFi BAM aligned to GRCh38.
#     pangu -p <prefix> --verbose [--vcf VCF] <inBam>
#   Output: <prefix>_report.json -> a JSON LIST; [0].diplotype == "CYP2D6 *4x2/*5".
#   We strip the leading "CYP2D6 " gene label and keep *N / xN tokens intact.
# GRCh38-only; needs the per-gene long-read BAM slice (3-sample sub-cohort only).
set -uo pipefail
sample=$1 gene=$2 build=$3 bundle=$4 out=$5
TOOL=pangu

emit() { # status diplotype [phenotype]
  printf 'tool\t%s\nsample\t%s\ngene\t%s\nbuild\t%s\nstatus\t%s\ndiplotype\t%s\nphenotype\t%s\n' \
    "$TOOL" "$sample" "$gene" "$build" "$1" "${2:-}" "${3:-}" > "$out"
}

# pangu is CYP2D6-only.
[ "$gene" = cyp2d6 ] || { emit NA ""; exit 0; }
# GRCh38-only (BAM must be aligned to GRCh38).
[ "$build" = grch38 ] || { emit NA ""; exit 0; }

bam="$bundle/${sample}.${gene}.longread.bam"
# Long-read slice only exists for the 3-sample sub-cohort; absent -> NA.
[ -s "$bam" ] || { emit NA ""; exit 0; }

work=$(mktemp -d)
prefix="$work/${sample}"
[ -s "${bam}.bai" ] || samtools index "$bam" >/dev/null 2>&1

pangu -p "$prefix" --verbose "$bam" >/dev/null 2>&1 \
  || { emit ERROR ""; rm -rf "$work"; exit 0; }

report="${prefix}_report.json"
[ -s "$report" ] || { emit ERROR ""; rm -rf "$work"; exit 0; }

diplo=$(python3 - "$report" <<'PY'
import sys, json, re
diplo = ""
try:
    d = json.load(open(sys.argv[1]))
    rec = d[0] if isinstance(d, list) and d else (d if isinstance(d, dict) else {})
    diplo = rec.get("diplotype", "") or ""
    # "CYP2D6 *4x2/*5" -> "*4x2/*5"; drop any leading gene label
    diplo = re.sub(r"^\s*CYP2D6\s+", "", diplo).strip()
    if diplo.lower() in ("none", "no_call", "."):
        diplo = ""
except Exception:
    pass
print(diplo.replace(" ", "_"))
PY
)
rm -rf "$work"

[ -n "$diplo" ] && emit OK "$diplo" "" || emit ERROR ""
