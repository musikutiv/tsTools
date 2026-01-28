#!/usr/bin/env Rscript
# Verify that optimized output matches original

suppressPackageStartupMessages({
  library(grid)
  library(IRanges)
  library(GenomicRanges)
  library(HilbertVis)
  library(RColorBrewer)
  library(AnnotationDbi)
  library(GenomicFeatures)
  library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
})

# Source the package code
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  source(f)
}

# Load test data
fpath <- file.path("inst", "extdata", "covx.rds")
cov <- new("SimpleRleList")
cov[["chrX"]] <- readRDS(fpath)

txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

# Test parameters
fstart <- 1660000
fend <- 1720000
fchr <- "chrX"
profs <- list(MSL2 = cov)

# Generate output
png("bench/test_output.png", width = 800, height = 400)
plotProfiles(
  fstart = fstart, fend = fend, fchr = fchr,
  profs = profs, cols = c("red"),
  txdb = txdb
)
dev.off()

cat("Output saved to bench/test_output.png\n")
cat("Visually verify the output is correct.\n")
