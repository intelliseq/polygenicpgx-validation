#!/usr/bin/env python3
"""
Genome-build handling test for polygenic (which is GRCh38-only for PGx).

For each (sample, gene) that has a GRCh37 Phase-3 slice, CrossMap-lift it to GRCh38
and re-run polygenic, then compare the lifted-GRCh37 call against the native-GRCh38
call. This quantifies whether the "lift then call" workaround preserves correctness,
and surfaces any strand/coordinate breakage (e.g. CYP2C19 rs3758581 / *38, where the
GRCh38 reference base itself defines an allele).

This is a standalone analysis (not part of the Nextflow per-tool fan-out). It writes
results/build_handling.tsv. Native-GRCh37 callers (Aldy/PyPGx/Cyrius/StellarPGx/
Stargazer) run Phase-3 directly in the main harness and are compared there.

Prereqs: CrossMap (pip install CrossMap), the hg19ToHg38 chain, a GRCh38 reference,
and the polygenic image / install on PATH (pgs-compute). Inputs are the work/data
slices produced by fetch_data.py for both builds.

  python3 benchmark/bin/liftover_grch37.py \
      --workdata benchmark/work/data --genes-yml benchmark/conf/genes.yml \
      --chain hg19ToHg38.over.chain.gz --ref GRCh38.fa --out benchmark/results/build_handling.tsv
"""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

import yaml

STAR = re.compile(r"\*\d+")


def sh(cmd, **kw):
    print("[lift]", " ".join(map(str, cmd)), file=sys.stderr)
    return subprocess.run(cmd, check=True, **kw)


def polygenic_call(vcf: Path, model: str) -> str:
    """Run pgs-compute and return the haplotype_id (or 'Indeterminate'/'ERROR')."""
    try:
        out = subprocess.run(
            ["pgstk", "pgs-compute", "--vcf", str(vcf), "--model", model, "--print"],
            capture_output=True, text=True, check=True).stdout
        h = json.loads(out).get("haplotype_model", {}).get("haplotypes", {})
        return h.get("haplotype_id") or h.get("call_filled") or "Indeterminate"
    except Exception:
        return "ERROR"


def lift(vcf37: Path, chain: str, ref: str, out_dir: Path) -> Path | None:
    """CrossMap a GRCh37 VCF to GRCh38; bgzip+index. Returns the lifted path or None."""
    lifted = out_dir / (vcf37.stem.replace(".grch37", ".lifted38") + ".vcf")
    try:
        sh(["CrossMap", "vcf", chain, str(vcf37), ref, str(lifted)])
        sh(["bgzip", "-f", str(lifted)])
        sh(["tabix", "-f", "-p", "vcf", str(lifted) + ".gz"])
        return Path(str(lifted) + ".gz")
    except subprocess.CalledProcessError:
        return None


def majors(call: str):
    return set(STAR.findall(call or ""))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdata", required=True, type=Path)
    ap.add_argument("--genes-yml", required=True, type=Path)
    ap.add_argument("--chain", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--models-dir", default="models/pgx")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    genes = yaml.safe_load(args.genes_yml.read_text())["genes"]
    rows = []
    # pair up native-grch38 and grch37 slices by (sample, gene)
    for vcf37 in sorted(args.workdata.glob("**/*-*.grch37.vcf.gz")):
        stem = vcf37.name.replace(".grch37.vcf.gz", "")
        sample, gene = stem.rsplit("-", 1)
        model = next(iter(Path(args.models_dir).glob(f"{gene}-pharmvar-*.yml")), None)
        if model is None:
            continue
        vcf38 = vcf37.with_name(f"{stem}.grch38.vcf.gz")
        native = polygenic_call(vcf38, str(model)) if vcf38.exists() else "MISSING"
        lifted_vcf = lift(vcf37, args.chain, args.ref, vcf37.parent)
        lifted = polygenic_call(lifted_vcf, str(model)) if lifted_vcf else "LIFT_FAIL"
        agree = int(bool(majors(native)) and majors(native) == majors(lifted))
        rows.append([sample, gene, native, lifted, str(agree)])

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        f.write("sample\tgene\tnative_grch38\tlifted_from_grch37\tmajor_alleles_agree\n")
        for r in sorted(rows):
            f.write("\t".join(r) + "\n")
    n_agree = sum(int(r[4]) for r in rows)
    print(f"[lift] {n_agree}/{len(rows)} (sample,gene) agree native-38 vs lifted-37 -> {args.out}")


if __name__ == "__main__":
    main()
