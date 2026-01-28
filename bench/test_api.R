#!/usr/bin/env Rscript
# Test script for plotProfiles API changes
# Verifies backwards compatibility and new GRanges API

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

# Test data setup
fpath <- file.path("inst", "extdata", "covx.rds")
cov <- new("SimpleRleList")
cov[["chrX"]] <- readRDS(fpath)
profs <- list(MSL2 = cov)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

# Test counter
tests_passed <- 0
tests_failed <- 0

test <- function(name, expr) {
  cat(sprintf("TEST: %s... ", name))
  result <- tryCatch({
    eval(expr)
    cat("PASS\n")
    tests_passed <<- tests_passed + 1
    TRUE
  }, error = function(e) {
    cat(sprintf("FAIL: %s\n", e$message))
    tests_failed <<- tests_failed + 1
    FALSE
  })
  invisible(result)
}

# Run tests
cat("\n========== API TESTS ==========\n\n")

# Test 1: Legacy API still works
test("Legacy API (fstart/fend/fchr)", {
  pdf(NULL); dev.control("enable")
  grob <- plotProfiles(
    fstart = 1660000, fend = 1720000, fchr = "chrX",
    profs = profs, cols = c("red"), txdb = txdb
  )
  dev.off()
  stopifnot(inherits(grob, "grob") || inherits(grob, "gTree"))
})

# Test 2: GRanges single range
test("GRanges API - single range", {
  pdf(NULL); dev.control("enable")
  ranges <- GRanges("chrX", IRanges(1660000, 1720000))
  grob <- plotProfiles(ranges = ranges, profs = profs, cols = c("red"), txdb = txdb)
  dev.off()
  stopifnot(inherits(grob, "grob") || inherits(grob, "gTree"))
})

# Test 3: GRanges multiple ranges
test("GRanges API - multiple ranges", {
  pdf(NULL); dev.control("enable")
  ranges <- GRanges("chrX", IRanges(c(1660000, 1720000), c(1720000, 1780000)))
  grobs <- plotProfiles(ranges = ranges, profs = profs, cols = c("red"), txdb = txdb)
  dev.off()
  stopifnot(is.list(grobs))
  stopifnot(length(grobs) == 2)
  stopifnot(inherits(grobs[[1]], "grob") || inherits(grobs[[1]], "gTree"))
})

# Test 4: Named ranges use names as titles
test("GRanges API - named ranges", {
  pdf(NULL); dev.control("enable")
  ranges <- GRanges("chrX", IRanges(c(1660000, 1720000), c(1720000, 1780000)))
  names(ranges) <- c("Region A", "Region B")
  grobs <- plotProfiles(ranges = ranges, profs = profs, cols = c("red"), txdb = txdb)
  dev.off()
  stopifnot(names(grobs)[1] == "Region A")
  stopifnot(names(grobs)[2] == "Region B")
})

# Test 5: Warning when both ranges and legacy args provided
test("Warning when ranges + legacy args provided", {
  pdf(NULL); dev.control("enable")
  ranges <- GRanges("chrX", IRanges(1660000, 1720000))
  w <- tryCatch({
    plotProfiles(
      ranges = ranges, fstart = 1660000, fend = 1720000, fchr = "chrX",
      profs = profs, cols = c("red"), txdb = txdb
    )
    NULL
  }, warning = function(w) w)
  dev.off()
  stopifnot(!is.null(w))
  stopifnot(grepl("ignoring deprecated", w$message))
})

# Test 6: Error when no range specified
test("Error when no range specified", {
  result <- tryCatch({
    plotProfiles(profs = profs, txdb = txdb)
    FALSE
  }, error = function(e) {
    grepl("Either 'ranges'", e$message)
  })
  stopifnot(result)
})

# Test 7: Error for invalid chromosome
test("Error for invalid chromosome", {
  pdf(NULL); dev.control("enable")
  result <- tryCatch({
    plotProfiles(fstart = 1000, fend = 2000, fchr = "chrZ",
                 profs = profs, cols = c("red"), txdb = txdb)
    FALSE
  }, error = function(e) {
    grepl("not found in coverage", e$message)
  })
  dev.off()
  stopifnot(result)
})

# Summary
cat(sprintf("\n========== RESULTS ==========\n"))
cat(sprintf("Passed: %d\n", tests_passed))
cat(sprintf("Failed: %d\n", tests_failed))
cat("=============================\n")

if (tests_failed > 0) {
  quit(status = 1)
}
