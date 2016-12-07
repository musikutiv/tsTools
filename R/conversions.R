#' Convert a bam file to coverage
#'
#' Reads can optionally be filtered by a chromosome inclusion list
#' Single reads will be extended to the provided fragment length
#'
#'
#' @param file.id path of the bam file
#' @param type sequencing library type "SINGLE" or "PAIRED"
#' @param fragment.length average chromatin fragment size, default 200
#' @param chr.flt vector of chromosome names to be included
#'
#'
#' @return a coverage object (\code{\link[IRanges]{RleList}} as returned by \code{\link[IRanges]{coverage}})
#'
#' @export
bam2coverage <- function(file.id, type=c("SINGLE","PAIRED"), fragment.length=200, chr.flt=NA) {

    if (!requireNamespace("ShortRead", quietly = TRUE)) {
    stop("library 'ShortRead' is needed for this function to work. Please install it.",
         call. = FALSE)
  }

    srbam2cov <- function(bam.file, extend=150, chr.flt=NA) {
    gr <- ShortRead::readGAlignments(bam.file)
    if (!is.na(chr.flt) & length(chr.flt)!=0) {
      gr <- keepSeqlevels(gr, chr.flt)
    }
    grs <- as(gr, "GRanges")
    grsr <- GRanges::resize(grs, extend)
    covs <- GRanges::coverage(grsr)
    return(covs)
  }

  # paired end to coverage
  # bam.file - BAM file
  # chr.flt - string vector of chromosome names to be included
  ## speed up!!!
  pebam2cov <- function(bam.file, chr.flt=NA) {
    gr <- readGAlignmentPairs(bam.file)
    flt <- seqnames(gr) %in% chr.flt
    gr <- gr[flt]
    gr <- grglist(gr)
    gr <- unlist(range(gr))
    covs <- coverage(gr)
    return(covs)
  }

  if (type=="SINGLE") {
    cov <- srbam2cov(file.id , chr.flt=chr.flt, extend=fragment.length)
  } else {
    cov <- pebam2cov(file.id , chr.flt=chr.flt)
  }
  cov
}



#' Convert numeric matrix to colors
#'
#' Useful for plotting heatmaps
#'
#' @param mat a matrix of values
#' @param rng mapping range
#'
#' @return a matrix of colors
#'
convertToColors <- function(mat, rng) {
  # Produce 'normalized' version of matrix, with values ranging from 0 to 1
  #rng <- range(mat, na.rm = TRUE)
  mat[mat<rng[1]] <- rng[1]
  mat[mat>rng[2]] <- rng[2]
  m <- (mat - rng[1])/diff(rng)
  m[is.nan(m)] <- NA
  # Convert to a matrix of sRGB color strings
  m2 <- m; class(m2) <- "character"
  #    m2[!is.na(m2)] <- rgb(colorRamp(rainbow(256))(m[!is.na(m)]), max = 255)
  m2[!is.na(m2)] <- rgb(colorRamp(c("white","midnightblue"))(m[!is.na(m)]), max = 255)
  m2[is.na(m2)] <- "transparent"
  return(m2)
}
