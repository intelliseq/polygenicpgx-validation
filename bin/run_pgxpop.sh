#!/usr/bin/env bash
# Runner contract (shared by every run_<tool>.sh):
#   args: <sample> <gene> <build> <bundle_dir> <out_raw>
#   writes <out_raw> as a 2-col TSV with keys: tool sample gene build status diplotype phenotype
#   status: OK (a call was produced) | NA (tool/build/gene unsupported or input missing) | ERROR
#
# PGxPOP: https://github.com/PharmGKB/PGxPOP   (image 1.0; defs frozen ~April 2020)
#   python bin/PGxPOP.py --vcf in.vcf.gz -g GENE [--phased] --build grch38 -o out.tsv
#   GRCh38-native (it bundles GRCh38 definitions). Output TSV columns include
#     sample_id  gene  diplotype  ...  phenotype  ...
#   Gene set: CFTR CYP2B6 CYP2C9 CYP2C19 CYP2D6 CYP3A5 CYP4F2 DPYD IFNL3 NUDT15
#             SLCO1B1 TPMT UGT1A1 VKORC1.
#   The slices are single-sample, phased for GRCh38; we pass --phased there.
set -uo pipefail
sample=$1 gene=$2 build=$3 bundle=$4 out=$5
TOOL=pgxpop

emit() { # status diplotype [phenotype]
  printf 'tool\t%s\nsample\t%s\ngene\t%s\nbuild\t%s\nstatus\t%s\ndiplotype\t%s\nphenotype\t%s\n' \
    "$TOOL" "$sample" "$gene" "$build" "$1" "${2:-}" "${3:-}" > "$out"
}

# Map our lowercase gene -> PGxPOP gene id. All ten harness genes are in PGxPOP's set.
declare -A GMAP=(
  [cyp2d6]=CYP2D6 [cyp2c19]=CYP2C19 [cyp2c9]=CYP2C9 [cyp3a5]=CYP3A5 [cyp2b6]=CYP2B6
  [slco1b1]=SLCO1B1 [vkorc1]=VKORC1 [tpmt]=TPMT [nudt15]=NUDT15 [dpyd]=DPYD
)
G=${GMAP[$gene]:-}
[ -n "$G" ] || { emit NA ""; exit 0; }

# PGxPOP bundles GRCh38 definitions; we run the GRCh38 path only (grch37 -> NA).
[ "$build" = grch38 ] || { emit NA ""; exit 0; }

vcf="$bundle/${sample}-${gene}.${build}.vcf.gz"
[ -s "$vcf" ] || { emit NA ""; exit 0; }

# PGxPOP wants a bgzipped + tabix-indexed VCF. The slices are single-sample; PGxPOP
# handles single-sample VCFs (it iterates samples). pytabix also wants `chr`-less or
# `chr`-prefixed contigs consistent with its defs; we pass the slice through unchanged
# (the bundle is GRCh38 chr-named to match PharmCAT/PGxPOP definitions).
work=$(mktemp -d)
in="$work/in.vcf.gz"
cp "$vcf" "$in"
if [ -s "${vcf}.tbi" ]; then cp "${vcf}.tbi" "$in.tbi"; else tabix -p vcf "$in" 2>/dev/null; fi

res="$work/pgxpop.tsv"
python3 "${PGXPOP_DIR:-/opt/PGxPOP}/bin/PGxPOP.py" \
  --vcf "$in" -g "$G" --phased --build grch38 -o "$res" \
  >/dev/null 2>&1 || { emit ERROR ""; rm -rf "$work"; exit 0; }

[ -s "$res" ] || { emit ERROR ""; rm -rf "$work"; exit 0; }

read -r diplo pheno < <(python3 - "$res" "$G" <<'PY'
import sys, csv
path, gene = sys.argv[1], sys.argv[2].upper()
diplo = pheno = ""
try:
    rows = list(csv.reader(open(path)))   # PGxPOP writes CSV despite the .tsv name
    hdr = [c.strip() for c in rows[0]]
    def col(*names):
        for n in names:
            for i, h in enumerate(hdr):
                if h.lower() == n.lower():
                    return i
        return None
    gi = col("gene")
    di = col("diplotype")
    pi = col("phenotype")
    # pick the row for the requested gene (single-sample slice -> usually one matching row)
    for r in rows[1:]:
        if gi is not None and gi < len(r) and r[gi].strip().upper() == gene:
            if di is not None and di < len(r): diplo = r[di].strip()
            if pi is not None and pi < len(r): pheno = r[pi].strip()
            break
except Exception:
    pass
if diplo.lower() in ("", "na", "nan", "none", "uncallable", "indeterminate"): diplo = ""
print(diplo.replace(" ", "_"), (pheno or "").replace(" ", "_"))
PY
)
rm -rf "$work"

[ -n "$diplo" ] && emit OK "$diplo" "$pheno" || emit ERROR ""
