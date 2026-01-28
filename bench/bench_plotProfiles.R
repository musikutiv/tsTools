#!/usr/bin/env Rscript
# Benchmark script for plotProfiles
# Run with: Rscript bench/bench_plotProfiles.R

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

# Set up reproducible test data
set.seed(42)

# Load the coverage data from the package
fpath <- file.path("inst", "extdata", "covx.rds")
cov <- new("SimpleRleList")
cov[["chrX"]] <- readRDS(fpath)

# Create annotation data
fpath_bed <- file.path("inst", "extdata", "msl2.bed")
peaks <- read.delim(fpath_bed, header = FALSE)
peaks[, 1] <- paste0("chr", peaks[, 1])
ann <- list("MSL2 peaks" = data.frame(
  chr = peaks[, 1],
  start = peaks[, 2],
  end = peaks[, 3],
  col = "blue"
))

# TxDb object
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

# Benchmark parameters
fstart <- 1660000
fend <- 1720000
fchr <- "chrX"

# Create GRanges objects for new API
single_range <- GRanges("chrX", IRanges(fstart, fend))
names(single_range) <- "Region 1"

# Multiple ranges
multi_range <- GRanges("chrX", IRanges(
  c(1660000, 1720000, 1780000),
  c(1720000, 1780000, 1840000)
))
names(multi_range) <- c("Region A", "Region B", "Region C")

# Create multiple profiles (typical use case)
profs_single <- list(MSL2 = cov)
profs_three <- list(MSL2 = cov, MSL3 = cov, MSL4 = cov)

# Warm-up run
cat("Warming up...\n")
pdf(NULL)
dev.control("enable")
plotProfiles(
  fstart = fstart, fend = fend, fchr = fchr,
  profs = profs_single, cols = c("red"),
  txdb = txdb
)
dev.off()

# Benchmark function
run_benchmark <- function(name, expr, n_runs = 5) {
  times <- numeric(n_runs)
  cat(sprintf("\nBenchmarking: %s (%d runs)\n", name, n_runs))

  for (i in seq_len(n_runs)) {
    pdf(NULL)
    dev.control("enable")
    start_time <- Sys.time()
    eval(expr)
    end_time <- Sys.time()
    dev.off()
    times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
    cat(sprintf("  Run %d: %.4f seconds\n", i, times[i]))
  }

  list(
    name = name,
    times = times,
    median = median(times),
    mean = mean(times),
    sd = sd(times),
    min = min(times),
    max = max(times)
  )
}

# Run benchmarks
results <- list()

# Legacy API tests
cat("\n========== LEGACY API (fstart/fend/fchr) ==========\n")

results$legacy_single <- run_benchmark(
  "Legacy: Single profile",
  quote(plotProfiles(
    fstart = fstart, fend = fend, fchr = fchr,
    profs = profs_single, cols = c("red"),
    txdb = txdb
  ))
)

results$legacy_three <- run_benchmark(
  "Legacy: Three profiles",
  quote(plotProfiles(
    fstart = fstart, fend = fend, fchr = fchr,
    profs = profs_three, cols = brewer.pal(3, "Set1"),
    txdb = txdb
  ))
)

# GRanges API tests
cat("\n========== NEW API (GRanges) ==========\n")

results$granges_single <- run_benchmark(
  "GRanges: Single range",
  quote(plotProfiles(
    ranges = single_range,
    profs = profs_single, cols = c("red"),
    txdb = txdb
  ))
)

results$granges_multi <- run_benchmark(
  "GRanges: Multi-range (3 ranges)",
  quote(plotProfiles(
    ranges = multi_range,
    profs = profs_single, cols = c("red"),
    txdb = txdb
  ))
)

# Calculate per-range time for multi-range
per_range_time <- results$granges_multi$median / 3

# Print summary
cat("\n\n========== BENCHMARK SUMMARY ==========\n")
cat("Legacy API:\n")
cat(sprintf("  %-35s: median=%.4fs\n", results$legacy_single$name, results$legacy_single$median))
cat(sprintf("  %-35s: median=%.4fs\n", results$legacy_three$name, results$legacy_three$median))
cat("\nGRanges API:\n")
cat(sprintf("  %-35s: median=%.4fs\n", results$granges_single$name, results$granges_single$median))
cat(sprintf("  %-35s: median=%.4fs (total), %.4fs (per range)\n",
            results$granges_multi$name, results$granges_multi$median, per_range_time))
cat("=========================================\n")

# Verify return types
cat("\n\n========== RETURN TYPE VERIFICATION ==========\n")
pdf(NULL)
dev.control("enable")

# Single range returns single grob
grob_single <- plotProfiles(
  ranges = single_range,
  profs = profs_single, cols = c("red"),
  txdb = txdb
)
cat(sprintf("Single range return type: %s\n", class(grob_single)[1]))

# Multiple ranges return list of grobs
grob_multi <- plotProfiles(
  ranges = multi_range,
  profs = profs_single, cols = c("red"),
  txdb = txdb
)
cat(sprintf("Multi-range return type: %s (length=%d)\n", class(grob_multi)[1], length(grob_multi)))
cat(sprintf("Multi-range names: %s\n", paste(names(grob_multi), collapse = ", ")))
cat(sprintf("Each element type: %s\n", class(grob_multi[[1]])[1]))

dev.off()
cat("=========================================\n")

# Save results
saveRDS(results, "bench/benchmark_results.rds")
cat("\nResults saved to bench/benchmark_results.rds\n")
