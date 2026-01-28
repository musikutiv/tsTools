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

# Create multiple profiles (typical use case)
profs_single <- list(MSL2 = cov)
profs_three <- list(MSL2 = cov, MSL3 = cov, MSL4 = cov)

# Warm-up run
cat("Warming up...\n")
pdf(NULL)  # Use null device for benchmarking
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

results$single_profile <- run_benchmark(
  "Single profile",
  quote(plotProfiles(
    fstart = fstart, fend = fend, fchr = fchr,
    profs = profs_single, cols = c("red"),
    txdb = txdb
  ))
)

results$three_profiles <- run_benchmark(
  "Three profiles",
  quote(plotProfiles(
    fstart = fstart, fend = fend, fchr = fchr,
    profs = profs_three, cols = brewer.pal(3, "Set1"),
    txdb = txdb
  ))
)

results$with_annotation <- run_benchmark(
  "Single profile + annotation",
  quote(plotProfiles(
    fstart = fstart, fend = fend, fchr = fchr,
    profs = profs_single, cols = c("red"),
    txdb = txdb, ann = ann
  ))
)

results$larger_window <- run_benchmark(
  "Larger window (120kb)",
  quote(plotProfiles(
    fstart = fstart, fend = fstart + 120000, fchr = fchr,
    profs = profs_single, cols = c("red"),
    txdb = txdb
  ))
)

# Print summary
cat("\n\n========== BENCHMARK SUMMARY ==========\n")
for (r in results) {
  cat(sprintf("%-30s: median=%.4fs, mean=%.4fs, sd=%.4fs\n",
              r$name, r$median, r$mean, r$sd))
}
cat("=========================================\n")

# Profiling with Rprof
cat("\n\n========== PROFILING (Single profile) ==========\n")
Rprof("bench/profile_output.out", memory.profiling = TRUE, line.profiling = TRUE)
pdf(NULL)
dev.control("enable")
for (i in 1:3) {
  plotProfiles(
    fstart = fstart, fend = fend, fchr = fchr,
    profs = profs_single, cols = c("red"),
    txdb = txdb
  )
}
dev.off()
Rprof(NULL)

# Summarize profile
cat("\nProfile summary (top 20 functions by total time):\n")
prof_summary <- summaryRprof("bench/profile_output.out")
print(head(prof_summary$by.total, 20))

cat("\n\nProfile summary (top 20 functions by self time):\n")
print(head(prof_summary$by.self, 20))

cat("\n\nTotal sampling time:", prof_summary$sampling.time, "seconds\n")

# Save results
saveRDS(list(results = results, profile = prof_summary),
        "bench/benchmark_results.rds")
cat("\nResults saved to bench/benchmark_results.rds\n")
