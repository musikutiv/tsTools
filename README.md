# tsTools

R utility functions for functional genomics data processing, primarily focused on ChIP-seq, MNase-seq, and related assays. The package provides coverage extraction, genomic visualization, nucleosome array analysis, and matrix manipulation utilities.

tsTools is a loose collection of helper functions. It is **not** a workflow system, an interactive genome browser, or a comprehensive analysis framework. It does not provide alignment, peak calling, or statistical testing functionality.

## Installation

Install from GitHub:

```r
# Install dependencies from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("HilbertVis", "IRanges", "GenomicRanges"))

# Install tsTools
install.packages("devtools")
devtools::install_github("musikutiv/tsTools")
```

### Dependencies

**Imports** (required):
- GenomicRanges, IRanges (Bioconductor)
- HilbertVis (Bioconductor)
- grid, grDevices, stats (base R)
- RColorBrewer

**Suggests** (optional, for specific functions):
- ShortRead: required for `bam2coverage()`
- data.table: required for `bed2dyad()`
- parallel: required for `*Parallel()` functions and `ocampo*()`
- zoo: required for `ocampo*()`
- TxDb.Dmelanogaster.UCSC.dm6.ensGene: example annotation database
- rtracklayer, lattice, testthat, knitr, rmarkdown

tsTools is not on CRAN or Bioconductor.

## Package Structure

```
R/
├── browser.R       # plotProfiles() - genomic coverage visualization
├── conversions.R   # bam2coverage() - BAM to coverage conversion
├── graphs.R        # cumPlot() - cumulative distribution plotting
├── nucleosomes.R   # Matrix utilities, bed2dyad(), plotRasterHeatmap()
├── ocampo.R        # ocampo(), ocampo2() - nucleosome array fitting
└── windows.R       # coverageWindows*() - coverage matrix extraction

bench/              # Benchmark scripts for plotProfiles performance testing
inst/extdata/       # Example data (Drosophila chrX coverage, BED files)
data/               # Package data (annotation example)
```

## Functional Overview

| Domain | Functions | Source |
|--------|-----------|--------|
| Genomic visualization | `plotProfiles` | browser.R |
| Coverage extraction | `coverageWindowsStranded`, `coverageWindowsCenteredStranded`, parallel variants | windows.R |
| Data conversion | `bam2coverage`, `bed2dyad` | conversions.R, nucleosomes.R |
| Nucleosome array analysis | `ocampo`, `ocampo2` | ocampo.R |
| Matrix utilities | `bin.matrix`, `bin.matrix.rows`, `norm.square`, `norm.sum`, `meanScale` | nucleosomes.R |
| Plotting helpers | `cumPlot`, `plotRasterHeatmap` | graphs.R, nucleosomes.R |

## Genomic Visualization

### plotProfiles()

Static, scriptable genome browser-style plots. Renders coverage profiles as filled polygons with gene model annotations from TxDb databases. Designed for publication figures, not interactive exploration.

**Requires**: A TxDb object for gene annotations. Coverage data as `SimpleRleList` objects (typically from `coverage()`).

**Recent changes**: Batch TxDb queries replaced per-gene database lookups, yielding ~4-5x performance improvement. GRanges input support added.

### plotProfiles() Reference

**Purpose**: Plot one or more genomic coverage profiles with gene annotations.

**Inputs**:

```r
plotProfiles(
  ranges = NULL,        # GRanges object (preferred)
  fstart = NULL,        # Deprecated: start position (bp)
  fend = NULL,          # Deprecated: end position (bp)
  fchr = NULL,          # Deprecated: chromosome name

  profs,                # Named list of SimpleRleList coverage objects
  cols = c(),           # Colors for profiles
  ann = NULL,           # Optional annotation list
  ylabel = "coverage",

  ylims = list(),       # Y-axis limits per profile

  txdb,                 # TxDb object for gene models
  ftitle = NA,
  collapse = TRUE,      # Collapse to longest transcript
  with.genes.highlited = c(),
  plot.labels = TRUE,
  grid = FALSE,
  with.average = FALSE
)
```

**GRanges input** (preferred):
```r
regions <- GRanges("chrX", IRanges(c(1660000, 1720000), c(1720000, 1780000)))
names(regions) <- c("Region A", "Region B")
plotProfiles(ranges = regions, profs = profs, txdb = txdb)
```

**Legacy input** (deprecated but supported):
```r
plotProfiles(fstart = 1660000, fend = 1720000, fchr = "chrX",
             profs = profs, txdb = txdb)
```

**Multi-range behavior**: When `ranges` contains multiple regions, one plot is drawn per range (each on a new page). Named ranges use names as plot titles.

**Return value**:
- Single range: `gTree` grob (invisible)
- Multiple ranges: Named list of `gTree` grobs

**Non-features**:
- No zooming or panning (static output only)
- No BAM/BigWig input (requires pre-computed coverage)
- No differential tracks or statistical overlays
- No PDF/PNG output handling (use standard R graphics devices)

## Coverage Extraction

Functions for extracting coverage values from `RleList` objects into matrices, with strand awareness.

### coverageWindowsStranded() / coverageWindowsStrandedParallel()

Extract coverage for arbitrary genomic windows. Minus-strand windows are automatically reversed.

```r
# windows: GRanges with strand information
# coverage: RleList from coverage()
mat <- coverageWindowsStranded(windows, coverage)
```

### coverageWindowsCenteredStranded() / coverageWindowsCenteredStrandedParallel()

Extract coverage centered on genomic positions (e.g., TSSs, binding sites).

```r
mat <- coverageWindowsCenteredStranded(centers, window.size = 1000, coverage)
```

Out-of-bounds windows are silently removed.

## Data Conversion

### bam2coverage()

Convert BAM files to coverage `RleList`. Requires ShortRead package.

```r
cov <- bam2coverage("file.bam", type = "PAIRED")
cov <- bam2coverage("file.bam", type = "SINGLE", fragment.length = 200)
```

### bed2dyad()

Convert BED files to nucleosome dyad coverage. For MNase-seq analysis. Requires data.table.

```r
dyad_cov <- bed2dyad("fragments.bed", type = "PAIRED")
```

## Nucleosome Array Analysis

### ocampo() / ocampo2()

Fit periodic nucleosome array patterns to coverage profiles. Based on correlation with Gaussian-modeled nucleosome arrays. Returns estimated spacing, phase shift, and correlation per region.

```r
# references: data.frame with chr, start, end, strand
results <- ocampo(coverage, references, beforeRef = 200, afterRef = 1000,
                  spacing.low = 130, spacing.high = 220)
```

`ocampo2()` uses cross-correlation (ccf) instead of explicit shift enumeration.

## Matrix Utilities

Row-wise normalization and binning functions for coverage matrices:

```r
norm.square(mat)      # L2 normalization per row
norm.sum(mat)         # Sum normalization per row
meanScale(mat)        # Mean centering per row
bin.matrix(mat, 50)   # Column binning, returns vector of means
bin.matrix.rows(mat, 50)  # Column binning, preserves rows
```

### plotRasterHeatmap()

Quick raster-based heatmap visualization of matrices. Uses 1st/99th percentile for color scaling.

```r
plotRasterHeatmap(mat)
```

## Examples

### Basic genome browser plot

```r
library(tsTools)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

# Load example coverage
cov <- new("SimpleRleList")
cov[["chrX"]] <- readRDS(system.file("extdata", "covx.rds", package = "tsTools"))

# Single region
plotProfiles(
  ranges = GRanges("chrX", IRanges(1660000, 1720000)),
  profs = list(MSL2 = cov),
  cols = "red",
  txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene
)
```

### Coverage matrix extraction

```r
# Define regions of interest
tss <- GRanges("chrX", IRanges(c(1000000, 2000000), width = 1),
               strand = c("+", "-"))
names(tss) <- c("gene1", "gene2")

# Extract coverage matrix centered on TSS
mat <- coverageWindowsCenteredStranded(tss, window.size = 2000, coverage = cov)
```

## Performance Considerations

`plotProfiles()` performance depends primarily on TxDb query overhead. The current implementation batches gene/transcript/exon queries per genomic window rather than per gene, reducing database round-trips.

For windows containing many genes (>50), expect sub-second plotting times on typical hardware. Very large windows or uncollapsed transcript models will increase rendering time.

The `*Parallel()` variants of coverage extraction functions provide modest speedups for large region sets when processing chromosomes in parallel.

## API Stability

tsTools is research software. The API may change between versions.

Current deprecations:
- `fstart`, `fend`, `fchr` arguments to `plotProfiles()` are deprecated in favor of `ranges`. They remain functional but emit no warning unless used alongside `ranges`.

## Limitations

- **plotProfiles()** requires pre-computed coverage (`SimpleRleList`). It does not read BAM/BigWig files directly.
- **Coverage extraction functions** silently drop out-of-bounds windows rather than erroring or padding.
- **bam2coverage()** loads entire BAM files into memory. Not suitable for very large files without chromosome filtering.
- **ocampo()** functions assume a specific nucleosome array model (Gaussian peaks at regular intervals). Results are correlation-based, not statistically rigorous.
- No unit tests are currently included in the package distribution.

## License

MIT (see DESCRIPTION)

## Citation

No formal citation exists for this package. If used in publications, cite the GitHub repository:

```
Straub, T. tsTools: Toolbox for Functional Genomics Data Processing in R.
https://github.com/musikutiv/tsTools
```
