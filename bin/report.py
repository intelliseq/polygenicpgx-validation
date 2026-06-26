#!/usr/bin/env python3
"""
Render the benchmark report (results/report.md) and the capability matrix
(results/capability_matrix.md) from concordance.tsv + the tool registry.

report.md:
  - per-tool x gene concordance (allele / diplotype %), with N
  - CYP2D6 structural stratum broken out (deletion / duplication / hybrid / tandem)
  - build-handling summary (grch38 vs grch37 rows)
  - coverage gaps (NA / ERROR / unmodelled)
capability_matrix.md:
  - input type, native GRCh37/38, CNV/hybrid, gene scope, group, license/notes
"""
import argparse
from collections import defaultdict
from pathlib import Path

import yaml


def pct(num, den):
    return f"{100*num/den:.0f}% ({num}/{den})" if den else "-"


def load_conc(p: Path):
    rows = []
    for line in p.read_text().splitlines()[1:]:
        c = line.split("\t")
        rows.append(dict(tool=c[0], sample=c[1], gene=c[2], build=c[3], status=c[4],
                         stratum=c[5], truth=c[6], call=c[7],
                         allele=int(c[8]), diplo=int(c[9])))
    return rows


def table(headers, rows):
    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join(["---"] * len(headers)) + "|"]
    out += ["| " + " | ".join(str(x) for x in r) + " |" for r in rows]
    return "\n".join(out)


def report_md(rows, genes):
    by_tool = sorted({r["tool"] for r in rows})
    gene_keys = list(genes.keys())
    g38 = [r for r in rows if r["build"] == "grch38"]

    # per-tool x gene diplotype concordance on GRCh38
    grid = []
    for t in by_tool:
        line = [t]
        for g in gene_keys:
            sub = [r for r in g38 if r["tool"] == t and r["gene"] == g]
            ok = [r for r in sub if r["status"] == "OK"]
            line.append(pct(sum(r["diplo"] for r in ok), len(sub)) if sub else "-")
        grid.append(line)

    # CYP2D6 structural strata
    strata = ["deletion", "homozygous_deletion", "duplication",
              "hybrid_68", "hybrid_13", "tandem_36", "tandem_36_hom", "adversarial"]
    cnv_rows = []
    for t in by_tool:
        line = [t]
        for st in strata:
            sub = [r for r in g38 if r["tool"] == t and r["gene"] == "cyp2d6" and r["stratum"] == st]
            line.append(pct(sum(r["diplo"] for r in sub), len(sub)) if sub else "-")
        cnv_rows.append(line)

    # build handling: per-tool grch38 vs grch37 diplotype rate
    build_rows = []
    for t in by_tool:
        line = [t]
        for b in ["grch38", "grch37"]:
            sub = [r for r in rows if r["tool"] == t and r["build"] == b]
            line.append(pct(sum(r["diplo"] for r in sub), len(sub)) if sub else "-")
        build_rows.append(line)

    # coverage gaps
    na = defaultdict(int)
    for r in rows:
        if r["status"] != "OK":
            na[(r["tool"], r["status"])] += 1

    # Headline: fair star-allele concordance over genes where major-allele matching is
    # meaningful — excludes CYP2D6 (structural, needs WGS), and VKORC1/DPYD (reported in
    # SNP-genotype / cDNA nomenclature the scorer can't reconcile with star alleles).
    CLEAN = ["cyp2c19", "cyp2c9", "cyp3a5", "cyp2b6", "slco1b1", "nudt15", "tpmt"]
    head = []
    for t in by_tool:
        sub = [r for r in g38 if r["tool"] == t and r["gene"] in CLEAN]
        ok = [r for r in sub if r["status"] == "OK"]
        if not ok:
            continue
        head.append([t, pct(sum(r["diplo"] for r in ok), len(sub)),
                     pct(sum(r["allele"] for r in ok), len(sub))])
    head.sort(key=lambda r: -float(r[1].split("%")[0]) if "%" in r[1] else 0)

    md = ["# PGx benchmark results\n",
          "_Concordance vs GeT-RM consensus (46 GeT-RM x 1000G samples). Diplotype-exact %, "
          "sub-allele-drift tolerant (major-allele match)._\n",
          "## Headline: star-allele concordance, GRCh38 (7 clean SNV genes)\n",
          "_Genes: CYP2C19/2C9/3A5/2B6, SLCO1B1, NUDT15, TPMT. Excludes CYP2D6 (structural) and "
          "VKORC1/DPYD (nomenclature) — see Caveats._\n",
          table(["tool", "diplotype %", "allele %"], head), "",
          "## Diplotype concordance by gene (GRCh38)\n",
          table(["tool"] + gene_keys, grid), "",
          "## CYP2D6 structural stratum (GRCh38)\n",
          table(["tool"] + strata, cnv_rows), "",
          "## Genome-build handling (diplotype %, all genes)\n",
          "_polygenic is GRCh38-only; its GRCh37 column reflects the liftover build-test "
          "(see liftover_grch37.py). Native-GRCh37 tools run Phase 3 directly._\n",
          table(["tool", "grch38", "grch37"], build_rows), "",
          "## Coverage gaps (non-OK calls)\n",
          table(["tool", "status", "count"],
                [[t, s, n] for (t, s), n in sorted(na.items())]) or "_none_", "",
          "## Caveats (read before interpreting)\n",
          "- **Gene-slice only** (by design): structural/CNV callers that need genome-wide depth "
          "or full HiFi WGS — Cyrius, StellarPGx, pangu, pb-StarPhase — return NA/ERROR on the "
          "gene-region slices. This is why CYP2D6 deletions/duplications/hybrids score ~0 for every "
          "tool: a slice cannot reveal a deleted allele (e.g. NA18861 `*29/*5` reads as `*29/*29`).\n",
          "- **VKORC1**: truth is the rs9923231 genotype (`GA/GG/AA`); tools emit `-1639G/-1639A` "
          "or rs-notation. Same variant, unreconciled nomenclature -> scores 0 spuriously.\n",
          "- **DPYD**: truth (Star Allele Search) uses `*1`/rs-tokens; some tools emit cDNA "
          "(`c.1601G>A`) -> nomenclature mismatch, not necessarily disagreement.\n",
          "- **T1K / SpecImmune** need gated PharmVar/IPD references (not baked into the images) -> NA. "
          "**deCYPher / Chinook** have no public CLI -> NA. **Stargazer** is GRCh37-only (contrast run).\n",
          "- **polygenic** has no VKORC1/TPMT model (cells shown as `-`) and is GRCh38-only; its "
          "GRCh37 column comes from the CrossMap liftover build-test (results/build_handling.tsv).\n"]
    return "\n".join(md)


def matrix_md(tools):
    headers = ["tool", "group", "input", "GRCh37", "GRCh38", "CNV/hybrid", "genes", "notes"]
    rows = []
    for name, t in tools.items():
        builds = t.get("builds", [])
        inp = t.get("input", "")
        inp = ",".join(inp) if isinstance(inp, list) else inp
        rows.append([
            name, t.get("group", ""), inp or "-",
            "✓" if "grch37" in builds else "✗",
            "✓" if "grch38" in builds else ("?" if t.get("group") == "C" else "✗"),
            "✓" if (t.get("cnv") or name in ("cyrius", "aldy", "stellarpgx", "pypgx",
                    "pbstarphase", "pangu", "decypher", "chinook", "specimmune", "stargazer")) else "✗",
            ",".join(t["genes"]) if t.get("genes") else "panel",
            t.get("note", t.get("allele_defs", "")),
        ])
    return "# PGx tool capability matrix\n\n" + table(headers, rows) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--concordance", required=True, type=Path)
    ap.add_argument("--calls", type=Path)
    ap.add_argument("--tool-versions", required=True, type=Path)
    ap.add_argument("--genes-yml", required=True, type=Path)
    ap.add_argument("--report", required=True, type=Path)
    ap.add_argument("--matrix", required=True, type=Path)
    args = ap.parse_args()

    genes = yaml.safe_load(args.genes_yml.read_text())["genes"]
    tools = yaml.safe_load(args.tool_versions.read_text())["tools"]
    rows = load_conc(args.concordance)

    args.report.write_text(report_md(rows, genes))
    args.matrix.write_text(matrix_md(tools))
    print(f"[report] {len(rows)} rows -> {args.report}, {args.matrix}")


if __name__ == "__main__":
    main()
