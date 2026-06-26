#!/usr/bin/env bash
# Runner contract (shared by every run_<tool>.sh):
#   args: <sample> <gene> <build> <bundle_dir> <out_raw>
#   writes <out_raw> as a 2-col TSV with keys: tool sample gene build status diplotype phenotype
#   status: OK (a call was produced) | NA (tool/build/gene unsupported or input missing) | ERROR
#
# PAnno: https://github.com/PreMedKB/PAnno   (pip install panno; image 0.2.0)
#   panno -s SAMPLE -i germline.vcf -p POP -o OUTDIR  ->  OUTDIR/<sample>.html
#   GRCh38 only. population (-p) is required; we default to EUR (override via PANNO_POP).
#   PAnno writes ONLY the HTML report (diplotypes are in-memory then rendered), so we
#   parse the requested gene's diplotype out of the rendered HTML diplotype table.
#   NOTE (verify): HTML scraping is best-effort and tied to PAnno's report layout; if the
#   layout changes, the parser may miss the gene and we report ERROR (ran, not parsed)
#   rather than a fake call.
set -uo pipefail
sample=$1 gene=$2 build=$3 bundle=$4 out=$5
TOOL=panno

emit() { # status diplotype [phenotype]
  printf 'tool\t%s\nsample\t%s\ngene\t%s\nbuild\t%s\nstatus\t%s\ndiplotype\t%s\nphenotype\t%s\n' \
    "$TOOL" "$sample" "$gene" "$build" "$1" "${2:-}" "${3:-}" > "$out"
}

# Map our lowercase gene -> PAnno gene symbol. PAnno covers 52 genes incl. all ten here.
declare -A GMAP=(
  [cyp2d6]=CYP2D6 [cyp2c19]=CYP2C19 [cyp2c9]=CYP2C9 [cyp3a5]=CYP3A5 [cyp2b6]=CYP2B6
  [slco1b1]=SLCO1B1 [vkorc1]=VKORC1 [tpmt]=TPMT [nudt15]=NUDT15 [dpyd]=DPYD
)
G=${GMAP[$gene]:-}
[ -n "$G" ] || { emit NA ""; exit 0; }

# PAnno's built-in definitions are GRCh38 only.
[ "$build" = grch38 ] || { emit NA ""; exit 0; }

vcf="$bundle/${sample}-${gene}.${build}.vcf.gz"
[ -s "$vcf" ] || { emit NA ""; exit 0; }

POP=${PANNO_POP:-EUR}

# PAnno reads the VCF directly; a gene-region slice is a valid germline VCF subset.
# (PAnno will simply find no variants for the other 51 genes; the requested gene is
#  covered by the slice.) Decompress in case PAnno's reader prefers plain VCF.
work=$(mktemp -d)
in="$work/in.vcf"
if ! bcftools view "$vcf" -O v -o "$in" 2>/dev/null; then
  zcat "$vcf" > "$in" 2>/dev/null || cp "$vcf" "$in"
fi

outdir="$work/report"
mkdir -p "$outdir"
panno -s "$sample" -i "$in" -p "$POP" -o "$outdir" >/dev/null 2>&1 \
  || { emit ERROR ""; rm -rf "$work"; exit 0; }

html=$(ls "$outdir/${sample}.html" "$outdir"/*.html 2>/dev/null | head -1)
[ -s "$html" ] || { emit ERROR ""; rm -rf "$work"; exit 0; }

# Scrape the gene's diplotype from the report. PAnno renders an inferred-diplotype table
# where the gene symbol is followed by its diplotype (e.g. "CYP2C19  *1/*2"). We pull the
# nearest *N.../rsN token pair on/after the gene cell and reconstruct "*1/*2".
read -r diplo pheno < <(python3 - "$html" "$G" <<'PY'
import sys, re, html as _html
path, gene = sys.argv[1], sys.argv[2].upper()
txt = open(path, encoding="utf-8", errors="ignore").read()
# strip tags to plain text, collapse whitespace
plain = _html.unescape(re.sub(r"(?is)<(script|style).*?</\1>", " ", txt))
plain = re.sub(r"(?s)<[^>]+>", " ", plain)
plain = re.sub(r"\s+", " ", plain)

diplo = pheno = ""
# Find the gene mention; capture a short window after it and grab a diplotype token.
# Diplotype forms: "*1/*2", "*1/*1", "rs...", "Reference/*X", "A/G"-style for SNP genes.
dip_re = re.compile(
    r"((?:\*\d+\w*|rs\d+|Reference|[ACGT]+)\s*[/|]\s*(?:\*\d+\w*|rs\d+|Reference|[ACGT]+))"
)
for m in re.finditer(re.escape(gene), plain):
    window = plain[m.end(): m.end() + 120]
    dm = dip_re.search(window)
    if dm:
        diplo = dm.group(1)
        diplo = re.sub(r"\s*([/|])\s*", r"/", diplo).strip()
        # phenotype: a metabolizer-style phrase nearby, best-effort
        pm = re.search(
            r"(Ultrarapid|Rapid|Normal|Intermediate|Poor|Likely\s+\w+|Indeterminate|"
            r"Increased|Decreased|Normal\s+Function|Possible\s+\w+)\s+(?:Metabolizer|Function|Risk)?",
            plain[m.end(): m.end() + 300], re.I)
        if pm:
            pheno = pm.group(0).strip()
        break

if diplo.lower() in ("", "na", "none", "indeterminate", "reference/reference"):
    # keep an explicit reference call only if it's clearly *1/*1-style; else blank
    if diplo.lower() in ("", "na", "none", "indeterminate"):
        diplo = ""
print(diplo.replace(" ", "_"), (pheno or "").replace(" ", "_"))
PY
)
rm -rf "$work"

[ -n "$diplo" ] && emit OK "$diplo" "$pheno" || emit ERROR ""
