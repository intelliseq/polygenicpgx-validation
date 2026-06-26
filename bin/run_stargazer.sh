#!/usr/bin/env bash
# Runner contract (shared by every run_<tool>.sh):
#   args: <sample> <gene> <build> <bundle_dir> <out_raw>
#   writes <out_raw> as a 2-col TSV with keys: tool sample gene build status diplotype phenotype
#   status: OK (a call was produced) | NA (tool/build/gene unsupported or input missing) | ERROR
#
# Stargazer: https://stargazer.gs.washington.edu/stargazerweb/
#   stargazer <dt> <gb> <tg> <vcf> <out> [--gdf FILE] [--cg STR]
#   Our harness uses Stargazer for GRCh37 ONLY (grch38 -> NA).
#   We run SNV genes from the VCF slice (wgs, VCF-only mode; no GDF/CNV).
#   Output: <out>/genotype.txt  ->  header: gene name status hap1_main hap2_main ... phenotype
#     diplotype = hap1_main "/" hap2_main when status == 'g' (genotyped).
# If the `stargazer` binary is not installed (academic-registration download), -> NA.
set -uo pipefail
sample=$1 gene=$2 build=$3 bundle=$4 out=$5
TOOL=stargazer

emit() { # status diplotype [phenotype]
  printf 'tool\t%s\nsample\t%s\ngene\t%s\nbuild\t%s\nstatus\t%s\ndiplotype\t%s\nphenotype\t%s\n' \
    "$TOOL" "$sample" "$gene" "$build" "$1" "${2:-}" "${3:-}" > "$out"
}

# Stargazer in this harness is GRCh37-only.
[ "$build" = grch37 ] || { emit NA ""; exit 0; }
GB=hg19

# Map our lowercase gene -> Stargazer target-gene id (lowercase symbol).
declare -A GMAP=(
  [cyp2d6]=cyp2d6 [cyp2c19]=cyp2c19 [cyp2c9]=cyp2c9 [cyp3a5]=cyp3a5 [cyp2b6]=cyp2b6
  [slco1b1]=slco1b1 [vkorc1]=vkorc1 [tpmt]=tpmt [nudt15]=nudt15 [dpyd]=dpyd
)
TG=${GMAP[$gene]:-}
[ -n "$TG" ] || { emit NA ""; exit 0; }

# Binary must be present (academic-gated download). Missing -> NA, not ERROR.
command -v stargazer >/dev/null 2>&1 || { emit NA ""; exit 0; }

vcf="$bundle/${sample}-${gene}.${build}.vcf.gz"
[ -s "$vcf" ] || { emit NA ""; exit 0; }

# Stargazer reads a plain (or bgzipped) VCF; an SV/GDF slice is optional. Decompress
# to be safe, and pass a control GDF only if the bundle happens to carry one.
work=$(mktemp -d)
invcf="$work/in.vcf"
if ! zcat "$vcf" > "$invcf" 2>/dev/null; then cp "$vcf" "$invcf"; fi

GDFARGS=()
gdf="$bundle/${sample}-${gene}.${build}.gdf"
if [ -s "$gdf" ]; then GDFARGS=(--gdf "$gdf" --cg vdr); fi

outdir="$work/proj"
stargazer wgs "$GB" "$TG" "$invcf" "$outdir" "${GDFARGS[@]}" >/dev/null 2>&1 \
  || { emit ERROR ""; rm -rf "$work"; exit 0; }

res="$outdir/genotype.txt"
[ -s "$res" ] || { emit ERROR ""; rm -rf "$work"; exit 0; }

read -r diplo pheno < <(python3 - "$res" <<'PY'
import sys, csv
diplo = pheno = ""
try:
    rows = list(csv.reader(open(sys.argv[1]), delimiter="\t"))
    hdr = [c.strip() for c in rows[0]]
    def idx(name):
        return next((i for i, h in enumerate(hdr) if h.lower() == name), None)
    si, h1, h2, pi = idx("status"), idx("hap1_main"), idx("hap2_main"), idx("phenotype")
    if rows[1:]:
        r = rows[1]
        get = lambda i: r[i].strip() if i is not None and i < len(r) else ""
        status = get(si)
        a, b = get(h1), get(h2)
        # genotyped row with real haplotypes -> "<a>/<b>"; '.' placeholders mean no call
        if a and b and a != "." and b != ".":
            diplo = f"{a}/{b}"
        pheno = get(pi)
        if pheno == ".":
            pheno = ""
except Exception:
    pass
print(diplo.replace(" ", "_"), pheno.replace(" ", "_"))
PY
)
rm -rf "$work"

[ -n "$diplo" ] && emit OK "$diplo" "$pheno" || emit ERROR ""
