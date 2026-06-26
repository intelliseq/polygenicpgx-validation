#!/usr/bin/env bash
# Runner contract (shared by every run_<tool>.sh):
#   args: <sample> <gene> <build> <bundle_dir> <out_raw>
#   writes <out_raw> as a 2-col TSV with keys: tool sample gene build status diplotype phenotype
#   status: OK (a call was produced) | NA (tool/build/gene unsupported or input missing) | ERROR
#
# PharmCAT pipeline: https://pharmcat.clinpgx.org/using/Running-PharmCAT-Pipeline/
#   `pharmcat_pipeline <vcf> -o <dir> -bf <base>` runs preprocessor + named-allele
#   matcher + phenotyper + reporter. We read the calls-only TSV (-reporterCallsOnlyTsv),
#   keyed by header (Gene / Source Diplotype / Phenotype), so column order is irrelevant.
#   GRCh38 only. PharmCAT does NOT call CYP2D6 from sequence (expects an outside call) -> NA.
set -uo pipefail
sample=$1 gene=$2 build=$3 bundle=$4 out=$5
TOOL=pharmcat

emit() { # status diplotype [phenotype]
  printf 'tool\t%s\nsample\t%s\ngene\t%s\nbuild\t%s\nstatus\t%s\ndiplotype\t%s\nphenotype\t%s\n' \
    "$TOOL" "$sample" "$gene" "$build" "$1" "${2:-}" "${3:-}" > "$out"
}

# PharmCAT is GRCh38-only.
[ "$build" = grch38 ] || { emit NA ""; exit 0; }

# PharmCAT does not call CYP2D6 from sequence (needs an outside call).
[ "$gene" = cyp2d6 ] && { emit NA ""; exit 0; }

# Map our lowercase gene -> PharmCAT gene id (uppercase canonical symbol).
declare -A GMAP=(
  [cyp2c19]=CYP2C19 [cyp2c9]=CYP2C9 [cyp3a5]=CYP3A5 [cyp2b6]=CYP2B6
  [slco1b1]=SLCO1B1 [vkorc1]=VKORC1 [tpmt]=TPMT [nudt15]=NUDT15 [dpyd]=DPYD
)
G=${GMAP[$gene]:-}
[ -n "$G" ] || { emit NA ""; exit 0; }   # gene not handled by PharmCAT here

vcf="$bundle/${sample}-${gene}.${build}.vcf.gz"
[ -s "$vcf" ] || { emit NA ""; exit 0; }

# Copy VCF locally (writable cwd) and run the pipeline (preprocessor + matcher +
# phenotyper). The earlier -bf/-reporterCallsOnlyTsv flags broke on 2.15.5; the plain
# run writes <base>.phenotype.json which carries the called diplotype + phenotype.
work=$(mktemp -d)
cp -f "$vcf" "$work/in.vcf.gz"
pharmcat_pipeline "$work/in.vcf.gz" -o "$work" >/dev/null 2>&1 \
  || { emit ERROR ""; rm -rf "$work"; exit 0; }

pjson=$(ls "$work"/*.phenotype.json 2>/dev/null | head -1)
[ -s "$pjson" ] || { emit ERROR ""; rm -rf "$work"; exit 0; }

read -r diplo pheno < <(python3 - "$pjson" "$G" <<'PY'
import json, sys
path, gene = sys.argv[1], sys.argv[2]
diplo = pheno = ""
try:
    d = json.load(open(path))
    for src in ("CPIC", "DPWG"):
        g = d.get("geneReports", {}).get(src, {}).get(gene)
        if not g:
            continue
        dips = g.get("sourceDiplotypes") or g.get("recommendationDiplotypes") or []
        if dips:
            diplo = dips[0].get("label", "") or ""
            ph = dips[0].get("phenotypes") or []
            pheno = ph[0] if ph else ""
            break
except Exception:
    pass
print((diplo or "").replace(" ", "_"), (pheno or "").replace(" ", "_"))
PY
)
rm -rf "$work"

[ -n "$diplo" ] && emit OK "$diplo" "$pheno" || emit NA ""
