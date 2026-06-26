#!/usr/bin/env python3
"""
Per-tool summary table for the README: how many panel genes each caller covers and its
star-allele concordance with the GeT-RM consensus on the calls it resolves (`matches/resolved`
+ %). GRCh38, resolved-call basis (a call is "resolved" when the tool emits a definite diplotype,
i.e. not Indeterminate / NA / ERROR). CYP2D6 is excluded (structural; not a fair SNV comparison),
and ursaPGx is omitted.

Usage: python3 make_summary.py --concordance results/concordance.tsv [--out results/summary.md]
"""
import argparse
from collections import defaultdict
from pathlib import Path

PANEL = ["cyp2c19", "cyp2c9", "cyp3a5", "cyp2b6", "slco1b1", "nudt15", "dpyd", "tpmt", "vkorc1"]
# VCF/GRCh38 callers compared head-to-head; ursaPGx omitted by request.
TOOLS = ["polygenic", "pharmcat", "pypgx", "panno", "pgxpop", "aldy"]
LABEL = {"polygenic": "polygenic", "pharmcat": "PharmCAT", "pypgx": "PyPGx",
         "panno": "PAnno", "pgxpop": "PGxPOP", "aldy": "Aldy"}


def load(p):
    rows = []
    for ln in Path(p).read_text().splitlines()[1:]:
        c = ln.split("\t")
        if c[3] == "grch38" and c[2] in PANEL:
            rows.append(dict(tool=c[0], gene=c[2], status=c[4],
                             call=c[7], diplo=int(c[9])))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--concordance", required=True)
    ap.add_argument("--out")
    args = ap.parse_args()
    rows = load(args.concordance)

    lines = ["| Tool | Genes covered | Concordance (resolved calls) |",
             "|------|:-------------:|:----------------------------:|"]
    table = []
    for t in TOOLS:
        sub = [r for r in rows if r["tool"] == t]
        if not sub:
            continue
        resolved = [r for r in sub if r["status"] == "OK" and r["call"]
                    and not r["call"].lower().startswith("indeterm")]
        genes = {r["gene"] for r in resolved}
        m = sum(r["diplo"] for r in resolved)
        n = len(resolved)
        pct = f"{100*m/n:.0f}% ({m}/{n})" if n else "-"
        table.append((t, len(genes), m, n, pct))
    # sort by concordance %, then coverage
    table.sort(key=lambda r: (-(r[2] / r[3] if r[3] else 0), -r[1]))
    for t, ng, m, n, pct in table:
        lines.append(f"| {LABEL[t]} | {ng}/{len(PANEL)} | {pct} |")

    md = "\n".join(lines) + "\n"
    print(md)
    if args.out:
        Path(args.out).write_text(md)
        print(f"[summary] wrote {args.out}")


if __name__ == "__main__":
    main()
