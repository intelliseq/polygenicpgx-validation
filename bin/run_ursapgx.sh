#!/usr/bin/env bash
# Runner contract (shared by every run_<tool>.sh):
#   args: <sample> <gene> <build> <bundle_dir> <out_raw>
#   writes <out_raw> as a 2-col TSV with keys: tool sample gene build status diplotype phenotype
#   status: OK (a call was produced) | NA (tool/build/gene unsupported or input missing) | ERROR
#
# ursaPGx: https://github.com/coriell-research/ursaPGx  (R/Bioconductor, image 1.0.0)
#   Needs a PHASED VCF; assigns diplotypes from PharmVar star-allele definitions.
#   Our bundle's gene slices are PHASED for GRCh38 only -> grch37 => NA.
#   We invoke a tool-specific R driver baked into the image (/opt/ursapgx_call.R):
#     Rscript ursapgx_call.R <phased_vcf_gz> <GENE_SYMBOL>  -> prints diplotype (e.g. *1|*2)
#   The driver prints "__UNSUPPORTED__" if ursaPGx doesn't model the gene -> NA.
#   ursaPGx has no phenotype layer (PharmVar definitions only) -> phenotype stays blank.
set -uo pipefail
sample=$1 gene=$2 build=$3 bundle=$4 out=$5
TOOL=ursapgx

emit() { # status diplotype [phenotype]
  printf 'tool\t%s\nsample\t%s\ngene\t%s\nbuild\t%s\nstatus\t%s\ndiplotype\t%s\nphenotype\t%s\n' \
    "$TOOL" "$sample" "$gene" "$build" "$1" "${2:-}" "${3:-}" > "$out"
}

# Map our lowercase gene -> ursaPGx gene symbol (uppercase canonical). CYP2D6 from a VCF
# is NOT supported by ursaPGx (it routes CYP2D6 through Cyrius/BAM), so we drop it -> NA.
declare -A GMAP=(
  [cyp2c19]=CYP2C19 [cyp2c9]=CYP2C9 [cyp3a5]=CYP3A5 [cyp2b6]=CYP2B6
  [slco1b1]=SLCO1B1 [vkorc1]=VKORC1 [tpmt]=TPMT [nudt15]=NUDT15 [dpyd]=DPYD
)
G=${GMAP[$gene]:-}
[ -n "$G" ] || { emit NA ""; exit 0; }

# GRCh38 only: the bundle's phased slices are GRCh38; ursaPGx needs phased input.
[ "$build" = grch38 ] || { emit NA ""; exit 0; }

vcf="$bundle/${sample}-${gene}.${build}.vcf.gz"
[ -s "$vcf" ] || { emit NA ""; exit 0; }

# ursaPGx (VariantAnnotation) wants a bgzipped + tabix-indexed VCF.
work=$(mktemp -d)
in="$work/in.vcf.gz"
cp "$vcf" "$in"
if [ -s "${vcf}.tbi" ]; then cp "${vcf}.tbi" "$in.tbi"; else tabix -p vcf "$in" 2>/dev/null; fi

driver=${URSAPGX_DRIVER:-/opt/ursapgx_call.R}
diplo=$(Rscript "$driver" "$in" "$G" 2>/dev/null)
rc=$?
rm -rf "$work"

if [ $rc -ne 0 ]; then emit ERROR ""; exit 0; fi

diplo=$(printf '%s' "$diplo" | tr -d '\r' | tail -1 | sed 's/[[:space:]]*$//')
case "$diplo" in
  __UNSUPPORTED__) emit NA ""; exit 0 ;;          # gene not modelled by ursaPGx
  "")              emit ERROR ""; exit 0 ;;        # ran, produced no call
  *)               emit OK "$diplo" "" ;;          # no phenotype layer in ursaPGx
esac
