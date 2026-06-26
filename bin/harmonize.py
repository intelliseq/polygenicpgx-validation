#!/usr/bin/env python3
"""
Aggregate per-tool .raw key/value files into one canonical calls.tsv.

Each run_<tool>.sh emits a 2-col TSV (keys: tool sample gene build status diplotype
phenotype). Harmonisation here is light: parse the raw diplotype into a canonical
allele representation that score.py can compare build- and tool-independently.

calls.tsv columns:
  tool sample gene build status diplotype_raw alleles cnv phenotype
    alleles : sorted comma-joined MAJOR star alleles (e.g. "*1,*4"); rs-tokens kept as-is
    cnv     : 1 if the call carries a structural/copy-number marker (xN, +, del/*5)
"""
import argparse
import re
from pathlib import Path

STAR = re.compile(r"\*\d+")          # major star allele
RS = re.compile(r"rs\d+")            # rs-token alleles (DPYD/VKORC1 style)
CN = re.compile(r"x\d+|\+|del", re.I)
VKORC = re.compile(r"^[AG]{2}$")     # VKORC1 reported as GG/GA/AA


def canonical(diplo: str):
    """Return (sorted major-allele list, cnv_flag)."""
    if not diplo:
        return [], 0
    d = diplo.replace("_", " ")
    if VKORC.match(d.strip()):
        return [d.strip()], 0
    alleles = STAR.findall(d) + RS.findall(d)
    cnv = 1 if CN.search(d) else 0
    return sorted(alleles), cnv


def read_raw(p: Path) -> dict:
    kv = {}
    for line in p.read_text().splitlines():
        if "\t" in line:
            k, v = line.split("\t", 1)
            kv[k] = v
    return kv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", nargs="+", required=True)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    rows = []
    for f in args.raw:
        kv = read_raw(Path(f))
        if not kv.get("tool"):
            continue
        alleles, cnv = canonical(kv.get("diplotype", ""))
        rows.append([
            kv.get("tool", ""), kv.get("sample", ""), kv.get("gene", ""),
            kv.get("build", ""), kv.get("status", ""), kv.get("diplotype", ""),
            ",".join(alleles), str(cnv), kv.get("phenotype", ""),
        ])

    with open(args.out, "w") as f:
        f.write("tool\tsample\tgene\tbuild\tstatus\tdiplotype_raw\talleles\tcnv\tphenotype\n")
        for r in sorted(rows):
            f.write("\t".join(r) + "\n")
    print(f"[harmonize] {len(rows)} calls -> {args.out}")


if __name__ == "__main__":
    main()
