# polygenicpgx-validation

Open, reproducible validation of the **polygenic** pharmacogenomic (PGx) star-allele caller against
established tools, on public reference samples with consensus ground-truth genotypes.

Everything needed to reproduce the comparison is here: a containerised Nextflow harness, one Docker
image per caller, the sample panel and the GeT-RM truth set. No result files are committed — run the
harness and regenerate them yourself.

## What is compared

46 Coriell reference samples that are both **GeT-RM–characterised** (CDC consensus genotypes) and in
the **1000 Genomes Project** (public WGS), across nine core CPIC pharmacogenes on **GRCh38**:
CYP2C19, CYP2C9, CYP3A5, CYP2B6, SLCO1B1, NUDT15, DPYD, TPMT, VKORC1. Each caller runs in a pinned
Docker image; native outputs are harmonised to a canonical diplotype and scored against the GeT-RM
consensus (major-allele match, sub-allele drift tolerated). Nomenclature is normalised before scoring
(e.g. `Reference` is treated as `*1`; VKORC1 `GA` / `C/T` / `-1639G/A` are compared as a variant-allele
dose) so that differences reflect genotype, not notation.

## Result

Concordance is reported on the **calls each tool resolves** — a call is *resolved* when the tool emits
a definite diplotype. polygenic deliberately returns **Indeterminate** instead of guessing when the
evidence is insufficient, so its denominator is smaller (it resolves fewer, lower-confidence calls)
while its resolved calls are more often correct. "Genes covered" is the number of the nine panel genes
for which the tool produced at least one call.

| Tool | Genes covered | Concordance (resolved calls) |
|------|:-------------:|:----------------------------:|
| polygenic | 9/9 | 95% (316/332) |
| PyPGx | 9/9 | 88% (342/390) |
| PAnno | 9/9 | 87% (340/390) |
| PharmCAT | 9/9 | 83% (322/390) |
| PGxPOP | 9/9 | 78% (303/390) |

_VCF input, GRCh38, resolved-call basis. Regenerate with `bin/make_summary.py` after a run._

### Callers evaluated on a different footing

Two of the established callers do not run on the same VCF/GRCh38 input, so they are reported separately
rather than mixed into the table above. Their figures reflect that footing and are not directly
comparable.

| Tool | Input / build | Genes covered | Concordance (resolved calls) |
|------|---------------|:-------------:|:----------------------------:|
| Aldy | BAM/CRAM, GRCh38 | 9/9 | 83% (299/359) |
| Stargazer | VCF, GRCh37 | 9/9 | 70% (233/335) |

Aldy is a BAM-based caller, run here on gene-region CRAM slices (its copy-number-neutral reference
region is included in each slice so it can normalise); its principal strength, structural CYP2D6 from
whole-genome read depth, is not exercised by this panel. Stargazer is GRCh37-only and was run on that
build. The harness also runs Cyrius and StellarPGx (whole-genome–depth CYP2D6 callers) and several
long-read callers; these require inputs outside this panel's scope.

A few properties of polygenic relevant to clinical use, all reproducible with this harness:

- **It does not make doubtful calls.** Where the evidence cannot resolve an allele, polygenic returns
  *Indeterminate* rather than a confident wrong diplotype.
- **Reference-defining alleles.** It correctly calls alleles defined by the GRCh38 reference base
  itself (e.g. CYP2C19 `*38`, rs3758581), which can otherwise silently mis-label carriers as
  wild-type.
- **Native GRCh38.** It operates on GRCh38 genomic +strand coordinates, avoiding the call-corrupting
  errors that can arise when GRCh37 data are lifted over.

## Reproduce

Requirements: Docker, Nextflow (≥22), `bcftools`/`samtools`/`tabix`, and a local GRCh38 reference
FASTA (for the BAM/CRAM callers). polygenic itself is the published PyPI package
[`polygenic-pgx`](https://pypi.org/project/polygenic-pgx/) (pinned in `docker/polygenic.Dockerfile`).

```bash
# 1. build the tool images (one per caller)
bash build_images.sh

# 2. point the harness at a GRCh38 reference (needed only for BAM/long-read callers)
export PGXBENCH_GRCH38_REF=/path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa

# 3. run (a single sample/gene first to smoke-test)
nextflow run main.nf -profile docker --samples NA12878 --genes cyp2c19 --tools polygenic,pharmcat
# ...or the full panel
nextflow run main.nf -profile docker

# 4. score + summarise
python3 bin/harmonize.py --raw results/raw/*.raw --out results/calls.tsv
python3 bin/score.py --calls results/calls.tsv --truth conf/truth.tsv --samples conf/samples.tsv --out results/concordance.tsv
python3 bin/make_summary.py --concordance results/concordance.tsv
```

`run_benchmark.sh` wraps the environment for convenience. The GRCh37→GRCh38 build test uses
`bin/liftover_grch37.py` (needs CrossMap and an `hg19ToHg38` chain).

## Layout

```
bin/        harness scripts: fetch_data, build_truth, harmonize, score, make_summary, run_<tool>.sh
conf/       genes.yml, samples.tsv, truth.tsv (cited), tool_versions.yml
docker/     one Dockerfile per caller
models/pgx/ polygenic's PGx models (PharmVar-derived; CPIC TPMT; VKORC1 rs9923231)
main.nf     nextflow.config   build_images.sh   run_benchmark.sh
```

## Data provenance

- **Samples / truth:** CDC GeT-RM consensus (Pratt et al., *J Mol Diagn* 2025, PMID 40122159; the
  source for VKORC1 + TPMT), the Star Allele Search annotation of the 1000G panel (Gharani et al.,
  *BMC Genomics* 2024, PMC10811916), and the Cyrius CYP2D6 truth set (Chen et al., *Pharmacogenomics J*
  2021, PMC7997805). `conf/truth.tsv` is small, derived, per-row-cited data; the raw source tables are
  not redistributed here (regenerate with `bin/build_truth.py` from the published files).
- **Genotypes:** 1000 Genomes 30× NYGC panel (GRCh38) and Phase 3 (GRCh37), fetched at runtime by
  `bin/fetch_data.py` (no genomic data is committed).
- **Allele definitions:** PharmVar; the TPMT and VKORC1 models use standard CPIC/dbSNP definitions and
  are not fitted to these samples.

## Competing interests

polygenic is developed by one of the authors (M.P.); M.P. and M.K. are founders of Intelliseq SA,
which develops it. This benchmark is published in full precisely so the comparison does not rest on
the authors' assertion: it uses external, independently published ground truth, public allele
definitions, default settings for every tool, and the same harness and scoring throughout, and can be
re-run by anyone.

## License

All rights reserved (Intelliseq).
