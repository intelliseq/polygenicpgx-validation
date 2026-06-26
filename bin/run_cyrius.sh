#!/usr/bin/env bash
# Runner contract (shared by every run_<tool>.sh):
#   args: <sample> <gene> <build> <bundle_dir> <out_raw>
#   writes <out_raw> as a 2-col TSV with keys: tool sample gene build status diplotype phenotype
#   status: OK (a call was produced) | NA (tool/build/gene unsupported or input missing) | ERROR
#
# Cyrius: https://github.com/Illumina/Cyrius  (v1.1.1)
#   CYP2D6 ONLY, from a WGS BAM/CRAM. Needs a manifest (one absolute CRAM path/line)
#   and the --genome build flag (19/37/38). For CRAM, pass --reference ref.fa.
#     python3 star_caller.py --manifest M --genome 38 --reference REF \
#             --prefix out --outDir D --threads 1
#   Output: <outDir>/out.tsv  ->  Sample <TAB> Genotype <TAB> Filter
#   Genotype e.g. "*1/*4", "*3/*68+*4", "*1x2/*2"; "None" = no-call.
# Cyrius really wants the full WGS coverage profile; we feed it the grch38 CRAM
# slice (the only BAM/CRAM the bundle carries). If no CRAM is present -> NA.
set -uo pipefail
sample=$1 gene=$2 build=$3 bundle=$4 out=$5
TOOL=cyrius

emit() { # status diplotype [phenotype]
  printf 'tool\t%s\nsample\t%s\ngene\t%s\nbuild\t%s\nstatus\t%s\ndiplotype\t%s\nphenotype\t%s\n' \
    "$TOOL" "$sample" "$gene" "$build" "$1" "${2:-}" "${3:-}" > "$out"
}

# Cyrius is CYP2D6-only.
[ "$gene" = cyp2d6 ] || { emit NA ""; exit 0; }

# Map our build label -> Cyrius --genome value.
case "$build" in
  grch38) GENOME=38 ;;
  grch37) GENOME=37 ;;
  *)      emit NA ""; exit 0 ;;
esac

# The bundle only carries a grch38 CRAM slice; Cyrius needs a BAM/CRAM (not VCF).
cram="$bundle/${sample}-${gene}.grch38.cram"
if [ "$build" != grch38 ] || [ ! -s "$cram" ]; then emit NA ""; exit 0; fi

# CRAM needs the reference FASTA (path provided by the harness env).
ref="${PGXBENCH_GRCH38_REF:-}"
[ -n "$ref" ] && [ -s "$ref" ] || { emit NA ""; exit 0; }

work=$(mktemp -d)
manifest="$work/manifest.txt"
printf '%s\n' "$cram" > "$manifest"

"${CYRIUS_DIR:-/opt/Cyrius}/star_caller.py" \
  --manifest "$manifest" --genome "$GENOME" --reference "$ref" \
  --prefix out --outDir "$work" --threads 1 >/dev/null 2>&1 \
  || python3 "${CYRIUS_DIR:-/opt/Cyrius}/star_caller.py" \
       --manifest "$manifest" --genome "$GENOME" --reference "$ref" \
       --prefix out --outDir "$work" --threads 1 >/dev/null 2>&1 \
  || { emit ERROR ""; rm -rf "$work"; exit 0; }

res="$work/out.tsv"
[ -s "$res" ] || { emit ERROR ""; rm -rf "$work"; exit 0; }

# Parse the single data row: column "Genotype". "None" / empty -> ERROR (ran, no call).
diplo=$(python3 - "$res" <<'PY'
import sys, csv
diplo = ""
try:
    rows = list(csv.reader(open(sys.argv[1]), delimiter="\t"))
    hdr = [c.strip() for c in rows[0]]
    gi = next((i for i, h in enumerate(hdr) if h.lower() == "genotype"), None)
    if gi is not None and rows[1:]:
        v = rows[1][gi].strip() if gi < len(rows[1]) else ""
        if v.lower() not in ("", "none", "no_call", "."):
            diplo = v
except Exception:
    pass
print(diplo.replace(" ", "_"))
PY
)
rm -rf "$work"

[ -n "$diplo" ] && emit OK "$diplo" "" || emit ERROR ""
