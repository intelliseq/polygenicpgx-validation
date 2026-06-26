#!/usr/bin/env Rscript
# ursaPGx star-allele driver for the PGx benchmark harness.
#   usage: Rscript ursapgx_call.R <phased_vcf_gz> <GENE_SYMBOL>
#   prints exactly one line to stdout: the diplotype (e.g. "*1|*2") or "" on no-call.
#   Any error -> non-zero exit (the bash runner treats that as ERROR).
#
# Docs: https://github.com/coriell-research/ursaPGx
# ursaPGx::callDiplotypes(vcf, gene=GENE, phased=TRUE) returns a DataFrame whose rows
# are samples and whose single column (named for the gene) holds the diplotype string.
# We feed a single-sample slice, so we take the first (only) row.
suppressMessages({
  library(ursaPGx)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("usage: ursapgx_call.R <phased_vcf_gz> <GENE>")
}
vcf  <- args[[1]]
gene <- toupper(args[[2]])

# Gate on the genes ursaPGx actually defines (PharmVar-backed). availableGenes() is the
# canonical list in the package; fall back to a known-supported set if the helper is
# absent in the installed version. CYP2D6 from VCF is not supported (needs Cyrius/BAM).
supported <- tryCatch(
  toupper(ursaPGx::availableGenes()),
  error = function(e) c(
    "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6", "CYP3A5",
    "CYP2B6", "CYP3A4", "CYP4F2", "SLCO1B1", "DPYD", "NUDT15", "TPMT", "VKORC1"
  )
)
if (!(gene %in% supported)) {
  # Not modelled by ursaPGx -> signal NA to the runner via a sentinel exit code.
  cat("__UNSUPPORTED__\n")
  quit(status = 0)
}

res <- tryCatch(
  ursaPGx::callDiplotypes(vcf, gene = gene, phased = TRUE),
  error = function(e) e
)
if (inherits(res, "error")) {
  message(conditionMessage(res))
  quit(status = 2)   # ran but failed -> ERROR
}

# Coerce the DataFrame/data.frame result to a plain matrix and pull the first cell.
df <- as.data.frame(res, stringsAsFactors = FALSE)
val <- ""
if (nrow(df) >= 1L && ncol(df) >= 1L) {
  # prefer the column named for the gene; otherwise first column
  col <- if (gene %in% toupper(colnames(df))) which(toupper(colnames(df)) == gene)[1] else 1L
  val <- as.character(df[1L, col])
}
if (is.na(val)) val <- ""
# normalise common no-call sentinels
if (val %in% c("NA", "None", "none", ".", "NA|NA", "NA/NA")) val <- ""
cat(val, "\n", sep = "")
