#' Convert a bed file to nucleosome dyad coverage
#'
#' Single reads will be shifted towards the dyad position assuming an end-localisation of the read
#' For paired reads the fragment mid points will be considered dyad position
#'
#' the bed files for SINGLE should have columns 1) chromosome 2) start 3) end 4) strand
#' the bed files for PAIRED should have columns 1) chromosome 2) start 3) end
#'
#' @param file.id path of the bed file
#' @param type sequencing library type "SINGLE" or "PAIRED"
#' @param width width of the dyad, default 1
#'
#' @return a coverage object (\code{\link[IRanges]{RleList}} as returned by \code{\link[IRanges]{coverage}})
#'
#' @export
bed2dyad <- function(file.id, type=c("SINGLE","PAIRED"), width=1) {

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("library 'data.table' is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  options(scipen=999)
  shift=73

  bed <- data.frame(data.table::fread(file.id, showProgress = T))

  if (type=="SINGLE") {
    bedR <- GenomicRanges::GRanges(bed[,1], IRanges::IRanges(bed[,2], bed[,3]), strand=bed[,4])
    dyads <- GenomicRanges::esize(bedR, width)
    amount <- strand(dyads)
    runValue(amount) <- ifelse(runValue(amount) == "-", -shift, shift)
    dyads <- GenomicRanges::shift(dyads, as.vector(amount))
  } else {
    bedR <- GenomicRanges::GRanges(bed[,1], IRanges::IRanges(bed[,2], bed[,3]))
    dyads <- GenomicRanges::resize(bedR, width, fix="center")
  }
  cov <- GenomicRanges::coverage(dyads)
  cov
}



#' Bin a matrix column-wise and average the bins across rows
#'
#' @param m a matrix
#' @param bin.size bin size
#'
#' @return a vector of binned values
#'
#' @export
bin.matrix <- function(m, bin.size) {
  bm <- c()
  bn <- c()
  for (i in 1:round((ncol(m)/bin.size))) {
    er <- (i*bin.size)
    sr <- er-(bin.size-1)
    bm <- append(bm,mean(m[,sr:er],na.rm=T))
    bn <- append(bn,round(mean(c(as.integer(colnames(m)[sr]),as.integer(colnames(m)[er])))))
  }
  names(bm) <- bn
  bm
}

#' Bin a matrix column-wise per row
#'
#' @param m a matrix
#' @param bin.size bin size
#'
#' @return a matrix of binned rowws
#'
#' @export
bin.matrix.rows <- function(m, bin.size) {
  splitv <- rep(1:round((ncol(m)/bin.size)), each=5)
  splitv <- append(splitv, rep(NA,ncol(m)-length(splitv)))
  bmat <- t(apply(m,1, function(x) {unlist(lapply(split(x,splitv), mean))}))
  cnames <- round(unlist(lapply(split(as.integer(colnames(m)), splitv),median)))
  colnames(bmat) <- cnames
}

#' Row-wise sum of squares normalization of a matrix
#'
#' @param mat a matrix
#'
#' @return a matrix of normalized values
#'
#' @export
norm.square <- function(mat) {t(apply(mat, 1, function(x){x/sqrt(sum(x^2))}))}


#' Plot Raster Heatmap of matrix
#'
#' @param mat the matrix to be visualized
#'
#'
#' @export
plotRasterHeatmap <- function(mat) {
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("library 'grid' is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  range <- c(quantile(mat, 0.01, na.rm=T),quantile(mat, 0.99, na.rm=T))
  grid::grid.newpage()
  grid::grid.raster(convertToColors(mat, range), height=unit(0.8, "npc"), width=unit(0.8, "npc"), interpolate = T)
}

