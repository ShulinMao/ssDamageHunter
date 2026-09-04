# Regenerates inst/extdata/example_candidates.csv.
#
# The shipped example is entirely synthetic: hand-chosen depths, alt-read counts
# and allele balances covering the three hypotheses the model discriminates.
# It contains no real cell, donor, or case identifier and no real sequencing
# data. Run from the repository root:
#
#     Rscript data-raw/make_example_candidates.R
#
row <- function(cell_id, muttype, mutsig, dp, scalt, ab, pass = FALSE)
  data.frame(cell_id, muttype, mutsig, dp, scalt, ab, pass, stringsAsFactors = FALSE)

d <- rbind(
  # --- single-stranded damage: af ~ 0.5 * ab -----------------------------------
  row("cell_A", "snv", "TCT:C>T", 60, 15, 0.50),
  row("cell_A", "snv", "ACA:C>A", 60, 14, 0.50),
  row("cell_A", "snv", "GTG:T>C", 80, 20, 0.50),
  row("cell_A", "snv", "CCG:C>G", 48, 12, 0.50),
  row("cell_A", "indel", "1:Del:C:1", 60, 16, 0.50),
  row("cell_B", "snv", "ATA:T>G", 72, 18, 0.50),
  row("cell_B", "snv", "GCC:C>T", 60, 13, 0.50),
  row("cell_B", "indel", "1:Ins:T:2", 56, 14, 0.50),
  row("cell_C", "snv", "TAT:A>T", 100, 25, 0.50),
  row("cell_C", "indel", "1:Del:A:1", 96, 24, 0.50),
  # allele-imbalanced loci: the expectation follows the closer haplotype
  row("cell_B", "snv", "TTA:T>A", 64, 21, 0.65),
  row("cell_C", "snv", "CAC:A>G", 128, 22, 0.65),

  # --- heterozygous germline-like: af ~ ab, the double-stranded model wins -----
  row("cell_A", "snv", "TCT:C>T", 60, 30, 0.50),
  row("cell_B", "snv", "ACC:C>T", 60, 31, 0.50),
  row("cell_C", "indel", "1:Del:T:2", 60, 29, 0.50),

  # --- low-frequency errors: af ~ 0.25 * ab, the error model wins --------------
  row("cell_A", "snv", "GGA:G>T", 60,  7, 0.50),
  row("cell_B", "snv", "CTC:T>C", 60,  8, 0.50),
  row("cell_C", "snv", "AAT:A>C", 80,  9, 0.50),
  row("cell_C", "indel", "1:Ins:C:1", 60, 7, 0.50),

  # --- already called as a somatic mutation by SCAN2 (pass = TRUE) -------------
  # Sits right at the single-stranded expectation but is excluded by !pass.
  row("cell_A", "snv", "TGT:G>A", 60, 15, 0.50, pass = TRUE),
  row("cell_C", "snv", "CGC:G>A", 72, 18, 0.50, pass = TRUE)
)

d$af  <- round(d$scalt / d$dp, 6)
d$chr <- "chr1"
d$pos <- seq(1000000L, by = 25000L, length.out = nrow(d))

# SNV mutsig is "<context>:<ref>><alt>", e.g. "TCT:C>T".
snv <- d$muttype == "snv"
d$refnt <- NA_character_
d$altnt <- NA_character_
d$refnt[snv] <- sub("^.*:", "", sub(">.*$", "", d$mutsig[snv]))
d$altnt[snv] <- sub("^.*>", "", d$mutsig[snv])

# Indel mutsig is the COSMIC ID channel; give each a matching ref/alt.
indel_alleles <- rbind(
  c("1:Del:C:1", "AC",  "A"),
  c("1:Ins:T:2", "G",   "GT"),
  c("1:Del:T:2", "CT",  "C"),
  c("1:Ins:C:1", "T",   "TC"),
  c("1:Del:A:1", "GA",  "G"))
for (i in seq_len(nrow(indel_alleles))) {
  hit <- !snv & d$mutsig == indel_alleles[i, 1]
  d$refnt[hit] <- indel_alleles[i, 2]
  d$altnt[hit] <- indel_alleles[i, 3]
}
stopifnot(!any(is.na(d$refnt)), !any(is.na(d$altnt)))

d <- d[, c("cell_id","chr","pos","refnt","altnt","muttype","mutsig",
           "dp","scalt","af","ab","pass")]
write.csv(d, "inst/extdata/example_candidates.csv", row.names = FALSE, quote = FALSE)
cat("wrote", nrow(d), "rows\n")
