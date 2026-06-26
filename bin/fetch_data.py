#!/usr/bin/env python3
"""
Fetch the input bundle for one (sample, gene, build) from public data.

Produces, under --out:
  <sample>-<gene>.<build>.vcf.gz (+.tbi)   gene-slice VCF (all VCF callers)
  <sample>-<gene>.<build>.cram   (+.crai)  gene-slice CRAM (BAM callers; GRCh38 only)
  <sample>.<gene>.longread.bam   (+.bai)    long-read slice (LR callers; if --longread yes)
  meta.txt                                  status flags consumed by the runners

Data sources:
  GRCh38 VCF : 1000G NYGC 30x phased SNV/INDEL/SV panel (carries symbolic <CN>/<DEL>
               records needed by polygenic CNV calling).
  GRCh37 VCF : 1000G Phase 3 release (low-coverage genotypes; no SV).
  CRAM       : 1000G 2504 high-coverage 30x CRAMs (ENA paths via sequence.index);
               needs a local GRCh38 reference (env PGXBENCH_GRCH38_REF) — skipped if unset.
  long-read  : GIAB (NA12878) / HPRC HiFi; resolved from env PGXBENCH_LONGREAD_<SAMPLE>.

Generalises scripts/build_pgx_fixtures.py. Remote slicing uses bcftools/samtools
(tabix range queries) so only the gene window is downloaded. Failures are non-fatal:
the bundle is still written with a status flag so downstream marks the cell ERROR/NA.

Prereqs on PATH: bcftools (>=1.16), tabix, samtools (for CRAM/BAM).
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

import yaml

PANEL_GRCH38 = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
    "1kGP_high_coverage_Illumina.chr{c}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
)
PANEL_GRCH37 = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
)
SEQ_INDEX = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index"
)


def run(cmd, **kw):
    print("[fetch]", " ".join(map(str, cmd)), file=sys.stderr)
    return subprocess.run(cmd, check=True, **kw)


def fetch_vcf(sample, gene, build, chrom, start, end, out: Path, cache: str) -> bool:
    url = (PANEL_GRCH38 if build == "grch38" else PANEL_GRCH37).format(c=chrom.replace("chr", ""))
    region = f"{chrom}:{start}-{end}"
    dest = out / f"{sample}-{gene}.{build}.vcf.gz"
    try:
        view = subprocess.Popen(
            ["bcftools", "view", "-Ou", "-s", sample, "-r", region, url],
            stdout=subprocess.PIPE, cwd=cache)
        strip = subprocess.Popen(
            ["bcftools", "annotate", "-x", "INFO/CSQ", "-Oz", "-o", str(dest.resolve())],
            stdin=view.stdout, cwd=cache)
        view.stdout.close()
        if strip.wait() != 0 or view.wait() != 0:
            return False
        run(["tabix", "-f", "-p", "vcf", str(dest)])
        return True
    except subprocess.CalledProcessError:
        return False


def cram_url_for(sample, cache) -> str | None:
    """Resolve the high-coverage CRAM URL for a sample from the sequence.index.

    The 1000G sequence.index is ~50 MB; cache it ONCE in a shared location
    (PGXBENCH_SEQ_INDEX, else <workdata>/sequence.index) instead of per-bundle.
    """
    idx = Path(os.environ.get(
        "PGXBENCH_SEQ_INDEX",
        Path(__file__).resolve().parent.parent / "work" / "sequence.index"))
    idx.parent.mkdir(parents=True, exist_ok=True)
    if not idx.exists() or idx.stat().st_size == 0:
        try:
            run(["bash", "-c", f"curl -sf {SEQ_INDEX} -o {idx}"])
        except subprocess.CalledProcessError:
            return None
    for line in idx.read_text(errors="ignore").splitlines():
        if line.startswith("#") or "\t" not in line:
            continue
        cols = line.split("\t")
        if (len(cols) > 9 and cols[9] == sample) or f"/{sample}.final.cram" in cols[0]:
            # index gives an ftp:// URL on ftp.sra.ebi.ac.uk; http is more reliable in htslib
            return cols[0].replace("ftp://", "http://")
    return None


def fetch_cram(sample, gene, chrom, start, end, out: Path, cache: str) -> bool:
    ref = os.environ.get("PGXBENCH_GRCH38_REF")
    if not ref or not Path(ref).exists():
        print("[fetch] no PGXBENCH_GRCH38_REF -> skip CRAM (BAM tools = NA)", file=sys.stderr)
        return False
    url = cram_url_for(sample, cache)
    if not url:
        return False
    dest = out / f"{sample}-{gene}.grch38.cram"
    # Aldy normalises coverage against a copy-number-neutral region (CYP2D8, chr22);
    # include it in the slice so BAM callers can run on genes off chr22.
    regions = [f"{chrom}:{start}-{end}"]
    if chrom.replace("chr", "") != "22":
        regions.append("chr22:42151472-42152258")
    try:
        run(["samtools", "view", "-T", ref, "-C", "-o", str(dest), url, *regions], cwd=cache)
        run(["samtools", "index", str(dest)])
        return True
    except subprocess.CalledProcessError:
        return False


def fetch_longread(sample, gene, chrom, start, end, out: Path) -> bool:
    url = os.environ.get(f"PGXBENCH_LONGREAD_{sample}")
    ref = os.environ.get("PGXBENCH_GRCH38_REF")
    if not url:
        print(f"[fetch] no PGXBENCH_LONGREAD_{sample} -> skip long-read", file=sys.stderr)
        return False
    dest = out / f"{sample}.{gene}.longread.bam"
    cmd = ["samtools", "view", "-b", "-o", str(dest), url, f"{chrom}:{start}-{end}"]
    if ref:
        cmd[2:2] = ["-T", ref]
    try:
        run(cmd)
        run(["samtools", "index", str(dest)])
        return True
    except subprocess.CalledProcessError:
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--build", required=True, choices=["grch38", "grch37"])
    ap.add_argument("--genes-yml", required=True)
    ap.add_argument("--longread", default="no")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    genes = yaml.safe_load(Path(args.genes_yml).read_text())["genes"]
    chrom, start, end = genes[args.gene][args.build]
    args.out.mkdir(parents=True, exist_ok=True)
    # bcftools/samtools write remote-index caches to cwd; keep them out of the bundle.
    cache_dir = args.out / ".cache"
    cache_dir.mkdir(exist_ok=True)
    cache = str(cache_dir)

    flags = {"vcf": False, "cram": False, "longread": False}
    flags["vcf"] = fetch_vcf(args.sample, args.gene, args.build, chrom, start, end, args.out, cache)
    if args.build == "grch38":
        flags["cram"] = fetch_cram(args.sample, args.gene, chrom, start, end, args.out, cache)
        if args.longread == "yes":
            flags["longread"] = fetch_longread(args.sample, args.gene, chrom, start, end, args.out)

    (args.out / "meta.txt").write_text(
        f"sample\t{args.sample}\ngene\t{args.gene}\nbuild\t{args.build}\n"
        + "".join(f"{k}\t{'OK' if v else 'MISSING'}\n" for k, v in flags.items()))
    print(f"[fetch] {args.sample}/{args.gene}/{args.build}: {flags}", file=sys.stderr)


if __name__ == "__main__":
    main()
