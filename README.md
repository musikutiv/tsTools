## Requirements

### Bioconductor packages HilbertVis, IRanges, GenomicRanges

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("HilbertVis", "IRanges", "GenomicRanges"))
```

#### Pandoc for building vignettes

```
remotes::install_github("cderv/pandoc")
library(pandoc)
# Install version
pandoc_install("2.7.3")
```

## Installation

```
install.packages("devtools")

library(devtools)
# if pandoc is present
install_github("musikutiv/tsTools", build_vignettes=T)
# otherwise
install_github("musikutiv/tsTools", build_vignettes=F)

```

To view the documentation (only if pandoc is installed) 

```
library(tsTools)
browseVignettes("tsTools”)
```
