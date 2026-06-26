#!/usr/bin/env python3
"""
Score harmonised calls against the GeT-RM consensus truth set.

Three concordance levels per call (substring/sub-allele tolerant, matching the
existing pgx_truth.py convention so PharmVar sub-allele drift is not penalised):
  allele_presence : every truth major-allele appears in the call
  diplotype_exact : the call's major-allele multiset equals truth's
  (phenotype is reported separately when tools emit it)

Output concordance.tsv is per-(tool,sample,gene,build) with match booleans + the
CYP2D6 SV stratum, so report.py can aggregate any way. Calls with status != OK are
recorded as not-concordant but kept (so coverage gaps are visible, not hidden).
"""
import argparse
import re
from collections import Counter
from pathlib import Path

STAR = re.compile(r"\*\d+")
RS = re.compile(r"rs\d+")
# wild-type expressed as "Reference"/"ref" (PharmCAT/PyPGx/PAnno DPYD etc.) == *1.
REFWORD = re.compile(r"\b(reference|ref|wild[- ]?type|wt)\b", re.I)


def major_multiset(diplo: str) -> Counter:
    # normalise nomenclature so "Reference" compares equal to "*1" (e.g. DPYD wild-type).
    d = REFWORD.sub("*1", (diplo or "").replace("_", " ")).strip()
    return Counter(STAR.findall(d) + RS.findall(d))


def vkorc1_dose(s: str):
    """Canonicalise any rs9923231 representation to the variant-allele dose (0/1/2).
    Handles GeT-RM 'GG/GA/AA', genomic 'C/T', '-1639G/-1639A', 'Reference/c.-1639G>A',
    'Reference/rs9923231', 'rs9923231 reference (C)/...'. Returns None if unparseable."""
    s = (s or "").strip().lower().replace("_", " ")
    if not s:
        return None
    if re.fullmatch(r"[ga]{2}", s):           # promoter genotype, G=ref A=variant
        return s.count("a")
    toks = re.split(r"[/|]", s)
    if len(toks) != 2:
        return None
    dose = 0
    for t in toks:
        t = t.strip()
        if t == "*h1":                        # Aldy VKORC1 haplotype *H1 carries rs9923231 (-1639A)
            dose += 1
            continue
        if re.fullmatch(r"\*h\d+", t):        # other Aldy VKORC1 haplotypes are reference
            continue
        if "reference" in t or t in ("g", "c", "-1639g") or t.endswith("(c)") or t.endswith("(g)"):
            continue                          # reference allele
        if ("variant" in t or "-1639a" in t or ">a" in t or t in ("a", "t")
                or t.endswith("(t)") or t.endswith("(a)") or t == "rs9923231"):
            dose += 1                         # variant allele
    return dose


def load_truth(p: Path) -> dict:
    out = {}
    for line in p.read_text().splitlines()[1:]:
        s, g, diplo, src = line.split("\t")
        out[(s, g)] = (diplo, src)
    return out


def load_samples(p: Path) -> dict:
    out = {}
    if not p or not Path(p).exists():
        return out
    for line in Path(p).read_text().splitlines()[1:]:
        c = line.split("\t")
        out[c[0]] = c[2]  # cyp2d6_sv_class
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True, type=Path)
    ap.add_argument("--truth", required=True, type=Path)
    ap.add_argument("--samples", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    truth = load_truth(args.truth)
    sv_class = load_samples(args.samples)

    out_rows = []
    for line in args.calls.read_text().splitlines()[1:]:
        tool, sample, gene, build, status, diplo_raw, alleles, cnv, pheno = line.split("\t")
        key = (sample, gene)
        if key not in truth:
            continue
        t_diplo, t_src = truth[key]
        ok = status == "OK"
        if gene == "vkorc1":                  # SNP genotype, compare variant-allele dose
            td, cd = vkorc1_dose(t_diplo), vkorc1_dose(diplo_raw)
            diplo_match = int(ok and td is not None and td == cd)
            allele_match = diplo_match
        else:
            tset = major_multiset(t_diplo)
            cset = major_multiset(diplo_raw)
            allele_match = int(ok and tset and set(tset) <= set(cset))
            diplo_match = int(ok and tset == cset and bool(tset))
        # CYP2D6 stratum: structural truth vs SNP-only
        stratum = sv_class.get(sample, "") if gene == "cyp2d6" else "snp"
        out_rows.append([
            tool, sample, gene, build, status, stratum,
            t_diplo, diplo_raw or "", str(allele_match), str(diplo_match),
        ])

    with open(args.out, "w") as f:
        f.write("tool\tsample\tgene\tbuild\tstatus\tstratum\t"
                "truth\tcall\tallele_match\tdiplotype_match\n")
        for r in sorted(out_rows):
            f.write("\t".join(r) + "\n")
    print(f"[score] {len(out_rows)} scored calls -> {args.out}")


if __name__ == "__main__":
    main()
