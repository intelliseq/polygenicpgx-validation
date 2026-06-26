#!/usr/bin/env bash
# Runner contract (shared by every run_<tool>.sh):
#   args: <sample> <gene> <build> <bundle_dir> <out_raw>
#   writes <out_raw> as a 2-col TSV with keys: tool sample gene build status diplotype phenotype
#   status: OK (a call was produced) | NA (tool/build/gene unsupported or input missing) | ERROR
#
# PyPGx: https://pypgx.readthedocs.io/en/latest/cli.html  (sbslee/pypgx)
#   `pypgx run-ngs-pipeline GENE <out_dir> --variants in.vcf.gz --assembly GRCh37|GRCh38`
#   -> <out_dir>/results.zip (SampleTable[Results]); data.tsv holds Genotype + Phenotype.
#   Supports GRCh37 and GRCh38. SV genes (e.g. CYP2D6) want a BAM; from a SNV-only VCF
#   slice PyPGx still returns an SNV-based diplotype, so we run the VCF path for all genes.
set -uo pipefail
sample=$1 gene=$2 build=$3 bundle=$4 out=$5
TOOL=pypgx

emit() { # status diplotype [phenotype]
  printf 'tool\t%s\nsample\t%s\ngene\t%s\nbuild\t%s\nstatus\t%s\ndiplotype\t%s\nphenotype\t%s\n' \
    "$TOOL" "$sample" "$gene" "$build" "$1" "${2:-}" "${3:-}" > "$out"
}

# Map our lowercase gene -> PyPGx gene id (uppercase canonical symbol).
declare -A GMAP=(
  [cyp2d6]=CYP2D6 [cyp2c19]=CYP2C19 [cyp2c9]=CYP2C9 [cyp3a5]=CYP3A5 [cyp2b6]=CYP2B6
  [slco1b1]=SLCO1B1 [vkorc1]=VKORC1 [tpmt]=TPMT [nudt15]=NUDT15 [dpyd]=DPYD
)
G=${GMAP[$gene]:-}
[ -n "$G" ] || { emit NA ""; exit 0; }

case "$build" in
  grch38) ASM=GRCh38 ;;
  grch37) ASM=GRCh37 ;;
  *)      emit NA ""; exit 0 ;;
esac

vcf="$bundle/${sample}-${gene}.${build}.vcf.gz"
[ -s "$vcf" ] || { emit NA ""; exit 0; }

# PyPGx requires a bgzipped + tabix-indexed VCF.
work=$(mktemp -d)
in="$work/in.vcf.gz"
cp "$vcf" "$in"
if [ -s "${vcf}.tbi" ]; then cp "${vcf}.tbi" "$in.tbi"; else tabix -p vcf "$in" 2>/dev/null; fi

outdir="$work/pipe"
pypgx run-ngs-pipeline "$G" "$outdir" --variants "$in" --assembly "$ASM" \
  >/dev/null 2>&1 || { emit ERROR ""; rm -rf "$work"; exit 0; }

res="$outdir/results.zip"
[ -s "$res" ] || { emit ERROR ""; rm -rf "$work"; exit 0; }

read -r diplo pheno < <(python3 - "$res" <<'PY'
import sys, zipfile, csv, io
diplo = pheno = ""
try:
    with zipfile.ZipFile(sys.argv[1]) as z:
        name = next(n for n in z.namelist() if n.endswith("data.tsv"))
        rows = list(csv.reader(io.StringIO(z.read(name).decode()), delimiter="\t"))
    hdr = [c.strip() for c in rows[0]]
    def col(*names):
        for n in names:
            for i, h in enumerate(hdr):
                if h.lower() == n.lower():
                    return i
        return None
    di = col("Genotype", "Diplotype")
    pi = col("Phenotype")
    if rows[1:]:
        r = rows[1]                      # single-sample slice -> one data row
        if di is not None and di < len(r): diplo = r[di].strip()
        if pi is not None and pi < len(r): pheno = r[pi].strip()
except Exception:
    pass
# normalise PyPGx "NA"/empty sentinels
if diplo.lower() in ("", "na", "nan", "indeterminate"): diplo = ""
print(diplo.replace(" ", "_"), (pheno or "").replace(" ", "_"))
PY
)
rm -rf "$work"

[ -n "$diplo" ] && emit OK "$diplo" "$pheno" || emit ERROR ""
