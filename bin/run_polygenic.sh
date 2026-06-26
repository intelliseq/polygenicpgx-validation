#!/usr/bin/env bash
# Runner contract (shared by every run_<tool>.sh):
#   args: <sample> <gene> <build> <bundle_dir> <out_raw>
#   writes <out_raw> as a 2-col TSV with keys: tool sample gene build diplotype phenotype status
#   status: OK (a call was produced) | NA (tool/build/gene unsupported or input missing) | ERROR
set -uo pipefail
sample=$1 gene=$2 build=$3 bundle=$4 out=$5
TOOL=polygenic

emit() { # status diplotype [phenotype]
  printf 'tool\t%s\nsample\t%s\ngene\t%s\nbuild\t%s\nstatus\t%s\ndiplotype\t%s\nphenotype\t%s\n' \
    "$TOOL" "$sample" "$gene" "$build" "$1" "${2:-}" "${3:-}" > "$out"
}

# polygenic PGx is GRCh38-only (GRCh37 is covered by the liftover build-test instead).
[ "$build" = grch38 ] || { emit NA ""; exit 0; }

src="$bundle/${sample}-${gene}.${build}.vcf.gz"
[ -s "$src" ] || { emit NA ""; exit 0; }
# Copy into the (writable) cwd so VcfAccessor can write its <vcf>.idx.db sidecar
# even when the staged input dir is read-only.
vcf="./$(basename "$src")"
cp -f "$src" "$vcf"; [ -s "$src.tbi" ] && cp -f "$src.tbi" "$vcf.tbi"

# match any model file for the gene (`<gene>-pharmvar-*.yml`, `<gene>-cpic-*.yml`, `<gene>-1.0.0.yml`)
model=$(ls "${POLYGENIC_MODELS_DIR:-models/pgx}/${gene}-"*.yml 2>/dev/null | head -1)
[ -n "$model" ] || { emit NA ""; exit 0; }   # unmodelled gene

# Write JSON to a file (the result can be large — passing it as an argv string
# overflows ARG_MAX, e.g. for CYP2D6). Parse the file, not an argument.
pgstk pgs-compute --vcf "$vcf" --model "$model" --print >./result.json 2>/dev/null \
    || { emit ERROR ""; exit 0; }

read -r diplo pheno < <(python3 - ./result.json <<'PY'
import json, sys
try:
    d = json.load(open(sys.argv[1]))
    h = d.get("haplotype_model", {}).get("haplotypes", {})
    dm = d.get("diplotype_model", {})            # single-SNP / compound genes (VKORC1, ...)
    call = h.get("haplotype_id") or h.get("call_filled") or h.get("match") or dm.get("diplotype") or ""
    pheno = (d.get("haplotype_model", {}).get("phenotype", "") or dm.get("category", "")
             or d.get("description", {}).get("phenotype", ""))
    print((call or "Indeterminate").replace(" ", "_"), (pheno or "").replace(" ", "_"))
except Exception:
    print(" ")
PY
)
[ -n "$diplo" ] && emit OK "$diplo" "$pheno" || emit ERROR ""
