#!/usr/bin/env bash
# Runner contract (shared by every run_<tool>.sh):
#   args: <sample> <gene> <build> <bundle_dir> <out_raw>
#   writes <out_raw> as a 2-col TSV with keys: tool sample gene build status diplotype phenotype
#   status: OK (a call was produced) | NA (tool/build/gene unsupported or input missing) | ERROR
#
# T1K: https://github.com/mourisl/T1K  (v1.0.5)
#   CYP2D6 only here, against a PharmVar allele-definition reference built into the
#   image ($T1K_REF_DIR/cyp2d6_dna_seq.fa). We extract reads from the grch38 CRAM
#   slice to FASTQ and run T1K with -1/-2 (avoids needing the genome-coord file):
#     run-t1k -1 r1.fq -2 r2.fq -f <ref>_dna_seq.fa --preset kir-wgs -o P --od D -t 1
#   Output: <P>_genotype.tsv  cols:
#     gene  num_alleles  allele1  abund1  qual1  allele2  abund2  qual2  [secondary]
#   (ignore alleles with quality <= 0). Diplotype = "<allele1>/<allele2>" (uses *N).
set -uo pipefail
sample=$1 gene=$2 build=$3 bundle=$4 out=$5
TOOL=t1k

emit() { # status diplotype [phenotype]
  printf 'tool\t%s\nsample\t%s\ngene\t%s\nbuild\t%s\nstatus\t%s\ndiplotype\t%s\nphenotype\t%s\n' \
    "$TOOL" "$sample" "$gene" "$build" "$1" "${2:-}" "${3:-}" > "$out"
}

# T1K wrapper is CYP2D6-only.
[ "$gene" = cyp2d6 ] || { emit NA ""; exit 0; }

# T1K consumes reads (BAM/CRAM/FASTQ); the bundle only carries a grch38 CRAM slice.
[ "$build" = grch38 ] || { emit NA ""; exit 0; }

cram="$bundle/${sample}-${gene}.grch38.cram"
[ -s "$cram" ] || { emit NA ""; exit 0; }

# PharmVar reference must have been built into the image (see t1k.Dockerfile).
refseq=$(ls "${T1K_REF_DIR:-/opt/t1kref}/${T1K_GENE:-cyp2d6}_dna_seq.fa" 2>/dev/null | head -1)
[ -n "$refseq" ] && [ -s "$refseq" ] || { emit NA ""; exit 0; }
refcoord=$(ls "${T1K_REF_DIR:-/opt/t1kref}/${T1K_GENE:-cyp2d6}_dna_coord.fa" 2>/dev/null | head -1)

ref="${PGXBENCH_GRCH38_REF:-}"
REFARG=()
[ -n "$ref" ] && [ -s "$ref" ] && REFARG=(--reference "$ref")

work=$(mktemp -d)
r1="$work/r1.fq"; r2="$work/r2.fq"; r0="$work/r0.fq"
# Extract paired reads from the slice; collate to keep mates together.
if samtools collate -u -O "${REFARG[@]}" "$cram" "$work/tmp.collate" 2>/dev/null \
     | samtools fastq "${REFARG[@]}" -1 "$r1" -2 "$r2" -0 /dev/null -s "$r0" -n - >/dev/null 2>&1; then :; fi

run_t1k() {
  CARG=(); [ -n "$refcoord" ] && CARG=(-c "$refcoord")
  if [ -s "$r1" ] && [ -s "$r2" ]; then
    run-t1k -1 "$r1" -2 "$r2" -f "$refseq" "${CARG[@]}" \
      --preset kir-wgs -o t1k --od "$work" -t 1 >/dev/null 2>&1
  elif [ -s "$r0" ]; then
    run-t1k -u "$r0" -f "$refseq" "${CARG[@]}" \
      --preset kir-wgs -o t1k --od "$work" -t 1 >/dev/null 2>&1
  else
    return 1
  fi
}

if ! run_t1k; then emit ERROR ""; rm -rf "$work"; exit 0; fi

res="$work/t1k_genotype.tsv"
[ -s "$res" ] || { emit ERROR ""; rm -rf "$work"; exit 0; }

diplo=$(python3 - "$res" <<'PY'
import sys
diplo = ""
try:
    rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1]) if l.strip()]
    # one line per gene: gene num_alleles allele1 abund1 qual1 allele2 abund2 qual2 ...
    row = rows[0] if rows else []
    def allele(name_i, qual_i):
        if name_i < len(row) and row[name_i].strip():
            try:
                q = float(row[qual_i]) if qual_i < len(row) and row[qual_i] != "" else 1
            except ValueError:
                q = 1
            a = row[name_i].strip()
            if a and a != "." and q > 0:
                # T1K reports "CYP2D6*4" or allele-series "CYP2D6*4ABC"; keep the *token.
                star = a.split("*", 1)
                return "*" + star[1] if len(star) == 2 else a
        return ""
    a1 = allele(2, 4)
    a2 = allele(5, 7)
    if a1 and a2:
        diplo = f"{a1}/{a2}"
    elif a1:
        diplo = f"{a1}/{a1}"   # homozygous (single reported allele)
except Exception:
    pass
print(diplo.replace(" ", "_"))
PY
)
rm -rf "$work"

[ -n "$diplo" ] && emit OK "$diplo" "" || emit ERROR ""
