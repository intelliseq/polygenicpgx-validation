#!/usr/bin/env bash
# Runner contract (shared by every run_<tool>.sh):
#   args: <sample> <gene> <build> <bundle_dir> <out_raw>
#   writes <out_raw> as a 2-col TSV with keys: tool sample gene build status diplotype phenotype
#   status: OK (a call was produced) | NA (tool/build/gene unsupported or input missing) | ERROR
#
# Aldy: https://github.com/0xTCG/aldy   https://aldy.readthedocs.io/
#   `aldy genotype -p illumina -g GENE --genome hg38|hg19 [-r ref.fa] -o out.aldy in`
#   Accepts BAM/CRAM/VCF (CRAM needs -r reference). The .aldy file carries a
#   `#Solution N: ...` header; the diplotype (major star alleles) is in the `Major`
#   column as e.g. "*1/*4". We prefer the CRAM slice (grch38) and fall back to VCF.
set -uo pipefail
sample=$1 gene=$2 build=$3 bundle=$4 out=$5
TOOL=aldy

emit() { # status diplotype [phenotype]
  printf 'tool\t%s\nsample\t%s\ngene\t%s\nbuild\t%s\nstatus\t%s\ndiplotype\t%s\nphenotype\t%s\n' \
    "$TOOL" "$sample" "$gene" "$build" "$1" "${2:-}" "${3:-}" > "$out"
}

# Aldy-supported genes (others -> NA). Map our lowercase id -> Aldy gene id.
declare -A GMAP=(
  [cyp2d6]=CYP2D6 [cyp2c19]=CYP2C19 [cyp2c9]=CYP2C9 [cyp3a5]=CYP3A5
  [dpyd]=DPYD [tpmt]=TPMT
  # (Aldy also does CYP2C8/CYP3A4/CYP4F2, not in our panel; cyp2b6/slco1b1/vkorc1/nudt15 unsupported)
)
G=${GMAP[$gene]:-}
[ -n "$G" ] || { emit NA ""; exit 0; }

# Map our build label -> Aldy --genome value.
case "$build" in
  grch38) GENOME=hg38 ;;
  grch37) GENOME=hg19 ;;
  *)      emit NA ""; exit 0 ;;
esac

cram="$bundle/${sample}-${gene}.grch38.cram"
vcf="$bundle/${sample}-${gene}.${build}.vcf.gz"

work=$(mktemp -d)
aldyout="$work/out.aldy"
rc=1
if [ "$build" = grch38 ] && [ -s "$cram" ]; then
  # CRAM requires a reference; REF_PATH/REF_FASTA may be provided by the harness env.
  REFARG=()
  [ -n "${ALDY_REFERENCE:-}" ] && REFARG=(-r "$ALDY_REFERENCE")
  aldy genotype -p illumina -g "$G" --genome "$GENOME" "${REFARG[@]}" \
    -o "$aldyout" "$cram" >/dev/null 2>&1 && rc=0
fi
if [ "$rc" -ne 0 ] && [ -s "$vcf" ]; then
  # VCF fallback (no reference needed for VCF input).
  aldy genotype -p illumina -g "$G" --genome "$GENOME" \
    -o "$aldyout" "$vcf" >/dev/null 2>&1 && rc=0
fi

if [ "$rc" -ne 0 ]; then
  # distinguish "no input present" (NA) from "tool failed on present input" (ERROR)
  if { [ "$build" = grch38 ] && [ -s "$cram" ]; } || [ -s "$vcf" ]; then
    emit ERROR ""
  else
    emit NA ""
  fi
  rm -rf "$work"; exit 0
fi

[ -s "$aldyout" ] || { emit ERROR ""; rm -rf "$work"; exit 0; }

# Parse the diplotype from the `#Solution N: <diplotype>` header line(s); take the
# first solution. Falls back to the Major column of the first data row.
diplo=$(python3 - "$aldyout" <<'PY'
import sys, re
diplo = ""
try:
    hdr = []
    data = []
    for line in open(sys.argv[1]):
        s = line.rstrip("\n")
        if not s:
            continue
        if s.startswith("#Solution") and not diplo:
            # e.g. "#Solution 1: *1/*4"  or  "#Solution 1: CYP2D6*1, CYP2D6*4"
            body = s.split(":", 1)[1].strip()
            diplo = body
            break
        if s.startswith("#Sample") or s.startswith("Sample"):
            hdr = s.lstrip("#").split("\t")
        elif not s.startswith("#"):
            data.append(s.split("\t"))
    if not diplo and hdr and data:
        try:
            mi = next(i for i, h in enumerate(hdr) if h.strip().lower() == "major")
            majors = sorted({r[mi].strip() for r in data if mi < len(r) and r[mi].strip()})
            if majors:
                diplo = "/".join(majors[:2]) if len(majors) > 1 else "/".join(majors * 2)
        except StopIteration:
            pass
except Exception:
    pass
# normalise "CYP2D6*1, CYP2D6*4" -> "*1/*4" but keep xN / + tandems intact
if diplo and "/" not in diplo and "," in diplo:
    diplo = "/".join(p.strip() for p in diplo.split(","))
diplo = re.sub(r"[A-Za-z0-9]+\*", "*", diplo)   # drop "CYP2D6" gene prefix before each *
print(diplo.replace(" ", "_"))
PY
)
rm -rf "$work"

[ -n "$diplo" ] && emit OK "$diplo" "" || emit ERROR ""
