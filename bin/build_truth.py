#!/usr/bin/env python3
"""
Build the benchmark truth set: conf/truth.tsv + conf/samples.tsv.

Truth provenance (the only sources of consensus genotypes — see benchmark/README.md
for the download URLs):
  - Star Allele Search Table S1 (Gharani 2024, PMC10811916): phased diplotypes for
    CYP2C19, CYP2C9, CYP3A5, CYP2B6, SLCO1B1, NUDT15, DPYD across the 3,202 NYGC-30x
    1000G samples. CYP2D6 in this table is mostly "Amb" (short reads cannot resolve
    CYP2D6 structural alleles), so CYP2D6 truth comes from Cyrius instead.
  - GeT-RM Consolidated PGx+HLA table (Pratt 2025, PMID 40122159): the ONLY source for
    VKORC1 (rs9923231 GG/GA/AA) and TPMT among GeT-RM samples.
  - Cyrius truth set (Chen 2021, PMC7997805): CYP2D6 diplotypes with structural variants
    resolved. The 40 structural shortlist calls are curated into CYP2D6_TRUTH below
    (verbatim from the Cyrius Table S1 / Star Allele Search cross-reference).

Raw source files live (gitignored) under benchmark/.sources/. Curated outputs
(conf/truth.tsv, conf/samples.tsv) ARE committed.

Run: python3 benchmark/bin/build_truth.py [--sources-dir benchmark/.sources]
"""
import argparse
import csv
from pathlib import Path

# ---- The 44-sample shortlist -------------------------------------------------
# CYP2D6 structural truth (Cyrius), with an SV class used to stratify scoring.
# Diplotypes use PharmVar nomenclature; "+" = tandem on one chromosome, "xN" = N copies.
CYP2D6_TRUTH = {
    # gene deletions (*5)
    "NA18861": ("*29/*5", "deletion"),
    "NA12873": ("*1/*5", "deletion"),
    "HG00276": ("*4/*5", "deletion"),
    "NA10831": ("*4/*5", "deletion"),
    "NA18855": ("*1/*5", "deletion"),
    "NA18868": ("*2/*5", "deletion"),
    "NA18945": ("*1/*5", "deletion"),
    "NA18992": ("*1/*5", "deletion"),
    "NA19035": ("*2/*5", "deletion"),
    "NA19317": ("*5/*5", "homozygous_deletion"),
    "HG03225": ("*56/*5", "deletion"),
    "HG03246": ("*43/*5", "deletion"),
    "HG03259": ("*106/*5", "deletion"),
    # duplications (xN)
    "NA19109": ("*2x2/*29", "duplication"),
    "NA19207": ("*10/*2x2", "duplication"),
    "NA19226": ("*2/*2x2", "duplication"),
    "NA19819": ("*2/*4x2", "duplication"),
    "NA19920": ("*1/*4x2", "duplication"),
    "HG00436": ("*2x2/*71", "duplication"),
    "HG00421": ("*10x2/*2", "duplication"),
    # *68 hybrids
    "NA12878": ("*3/*68+*4", "hybrid_68"),
    "NA11832": ("*1/*68+*4", "hybrid_68"),
    "NA12154": ("*33/*68+*4", "hybrid_68"),
    "HG00731": ("*4/*68+*4", "hybrid_68"),
    "HG01190": ("*5/*68+*4", "hybrid_68_deletion"),
    # *36 tandems
    "NA18526": ("*1/*36+*36+*10", "tandem_36"),
    "NA18565": ("*36/*36+*10", "tandem_36"),
    "NA18545": ("*36+*10/*36+*10", "tandem_36_hom"),
    "NA18563": ("*1/*36+*10", "tandem_36"),
    "NA18564": ("*2/*36+*10", "tandem_36"),
    "NA18572": ("*41/*36+*10", "tandem_36"),
    "NA18617": ("*36+*10/*36+*10", "tandem_36_hom"),
    "NA18632": ("*52/*36+*36+*10", "tandem_36"),
    "NA18642": ("*1+*90/*36+*10", "tandem_36"),
    "NA18959": ("*2/*36+*10", "tandem_36"),
    "NA18980": ("*2/*36+*10", "tandem_36"),
    "HG00463": ("*36+*10/*36+*10", "tandem_36_hom"),
    "HG02373": ("*14/*36+*10", "tandem_36"),
    # *13 hybrid + adversarial
    "NA19785": ("*13+*2/*1", "hybrid_13"),
    "NA18519": ("*106/*29", "adversarial"),  # Cyrius corrected GeT-RM *1 -> *106
}

# Population labels (1000G), plus 6 SNP-reference samples that complete the other genes.
POPULATION = {
    "NA12878": "CEU", "NA18861": "YRI", "NA12873": "CEU", "HG00276": "FIN",
    "NA10831": "CEU", "NA18855": "YRI", "NA18868": "YRI", "NA18945": "JPT",
    "NA18992": "JPT", "NA19035": "YRI", "NA19317": "LWK", "HG03225": "ESN",
    "HG03246": "ESN", "HG03259": "ESN", "NA19109": "YRI", "NA19207": "ASW",
    "NA19226": "ASW", "NA19819": "ASW", "NA19920": "ASW", "HG00436": "CHS",
    "HG00421": "CHS", "NA11832": "CEU", "NA12154": "CEU", "HG00731": "PUR",
    "HG01190": "PUR", "NA18526": "CHB", "NA18565": "CHB", "NA18545": "CHB",
    "NA18563": "CHB", "NA18564": "CHB", "NA18572": "CHB", "NA18617": "CHB",
    "NA18632": "CHB", "NA18642": "CHB", "NA18959": "JPT", "NA18980": "JPT",
    "HG00463": "CHB", "HG02373": "ACB", "NA19785": "MXL", "NA18519": "YRI",
    # SNP-reference extras
    "NA19143": "ASW", "NA07000": "CEU", "NA19122": "YRI", "NA19238": "YRI",
    "NA19239": "YRI", "NA12892": "CEU",
}

# Long-read sub-cohort (3 CYP2D6 samples with public PacBio HiFi / ONT data).
# NA12878 = GIAB; the two HPRC *68-hybrid trio members confirmed at fetch time.
LONGREAD_COHORT = {"NA12878", "HG00731", "HG01190"}

# Genes whose truth comes from Star Allele Search Table S1 (column name -> gene key).
SAS_GENES = {
    "CYP2C19": "cyp2c19", "CYP2C9": "cyp2c9", "CYP3A5": "cyp3a5",
    "CYP2B6": "cyp2b6", "SLCO1B1": "slco1b1", "NUDT15": "nudt15", "DPYD": "dpyd",
}
# VKORC1 + TPMT only exist in the GeT-RM consolidated dump.
GETRM_GENES = {"VKORC1": "vkorc1", "TPMT": "tpmt"}


def _norm(diplo: str) -> str:
    """Normalise a source diplotype cell; '.'/''/'Amb' -> '' (no truth)."""
    d = (diplo or "").strip()
    if d in (".", "", "Amb", "*Amb|*Amb", "n/a", "NA"):
        return ""
    return d.replace("|", "/")  # phased a|b -> a/b for build-agnostic comparison


def load_sas(path: Path) -> dict:
    """sample -> {gene_key: diplotype} from Star Allele Search Table S1."""
    out = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            s = row["Sample"].strip()
            out[s] = {g: _norm(row.get(col, "")) for col, g in SAS_GENES.items()}
    return out


def load_getrm(path: Path) -> dict:
    """sample -> {vkorc1, tpmt} from the GeT-RM consolidated dump (TSV)."""
    out = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            s = row["Sample"].strip()
            out[s] = {g: _norm(row.get(col, "")) for col, g in GETRM_GENES.items()}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sources-dir", default="benchmark/.sources", type=Path)
    ap.add_argument("--out-dir", default="benchmark/conf", type=Path)
    args = ap.parse_args()

    sas = load_sas(args.sources_dir / "star_allele_search_TableS1.csv")
    getrm = load_getrm(args.sources_dir / "getrm_consolidated_dump.tsv")

    samples = sorted(POPULATION, key=lambda s: (s[:2], s))
    truth_rows = []
    for s in samples:
        # CYP2D6 (Cyrius, authoritative for structural)
        if s in CYP2D6_TRUTH:
            diplo, _cls = CYP2D6_TRUTH[s]
            truth_rows.append((s, "cyp2d6", diplo, "Cyrius_2021_PMC7997805"))
        # 7 SAS genes
        for gene, diplo in sas.get(s, {}).items():
            if diplo:
                truth_rows.append((s, gene, diplo, "StarAlleleSearch_2024_PMC10811916"))
        # VKORC1 + TPMT
        for gene, diplo in getrm.get(s, {}).items():
            if diplo:
                truth_rows.append((s, gene, diplo, "GeTRM_Consolidated_2025_PMID40122159"))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    truth_path = args.out_dir / "truth.tsv"
    with open(truth_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "gene", "diplotype", "source"])
        for r in sorted(truth_rows):
            w.writerow(r)

    samples_path = args.out_dir / "samples.tsv"
    with open(samples_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "population", "cyp2d6_sv_class", "longread"])
        for s in samples:
            cls = CYP2D6_TRUTH.get(s, ("", "snp_reference"))[1]
            w.writerow([s, POPULATION[s], cls, "yes" if s in LONGREAD_COHORT else "no"])

    print(f"[truth] {len(truth_rows)} rows across {len(samples)} samples -> {truth_path}")
    print(f"[samples] {len(samples)} samples -> {samples_path}")


if __name__ == "__main__":
    main()
