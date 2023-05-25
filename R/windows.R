
#' Get windowed strand-oriented coverages around center points
#'
#' Out-of-bound windows will be removed!
#'
#' @param centers a named GRanges object (\code{\link[GenomicRanges]{Granges}}
#' @param window.size the size of the window surrounding the center position, default 1000
#' @param coverage a coverage object (\code{\link[IRanges]{RleList}} as returned by \code{\link[IRanges]{coverage}})
#'
#' @return a matrix
#'
#' @export
coverageWindowsCenteredStranded <- function (centers, window.size = 1000, coverage) {

  ## only coverages with enters in centers
  coverage <- coverage[names(coverage) %in% seqlevels(centers)]

  # centers <- centers[centers$chr %in% names(coverage),]

  result <- lapply(names(coverage), function(x) {
    my.cov <- coverage[[x]]
    my.centers <- centers[seqnames(centers) == x, ]
    mw.views <- IRanges::Views(my.cov, start = start(my.centers) -
                                 ceiling(window.size/2), width = window.size + 1)
    flt <- start(mw.views) > 0 & end(mw.views) < length(my.cov)
    mw.views <- mw.views[flt, ]
    my.centers <- my.centers[flt]
    if (length(mw.views) > 0) {
      mat <- as.matrix(mw.views)
      colnames(mat) <- seq(from = (0 - ceiling(window.size/2)),
                           to = 0 + ceiling(window.size/2))
      rownames(mat) <- names(my.centers)
      return(mat)
    }
    else {
      return(NULL)
    }
  })
  mat <- Reduce(rbind, result)
  centers <- centers[names(centers) %in% rownames(mat),]
  o <- match(names(centers), rownames(mat))
  mat <- mat[o, ]
  if (sum(strand(centers) == "-") > 0) {
    mat[as.vector(strand(centers) == "-"), ] <- t(apply(mat[as.vector(strand(centers) == "-"), ], 1, rev))
  }
  colnames(mat) <- seq(-ceiling(window.size/2), ceiling(window.size/2))
  mat
}

#' Get windowed strand-oriented coverages around center point, parallelized
#'
#' Out-of-bound windows will be removed!
#'
#' @param centers a named GRanges object (\code{\link[GenomicRanges]{Granges}}
#' @param window.size the size of the window surrounding the center position, default 1000
#' @param coverage a coverage object (\code{\link[IRanges]{RleList}} as returned by \code{\link[IRanges]{coverage}})
#' @param n.cores number of processor cores to use. If not defined all detected cores will be used)
#'
#' @return a matrix
#'
#' @export
coverageWindowsCenteredStrandedParallel <- function(centers, window.size=1000, coverage, n.cores=2) {

  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("library 'parallel' is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (is.null(n.cores)) n.cores <- detectCores()

  ## only coverages with enters in centers
  coverage <- coverage[names(coverage) %in% seqlevels(centers)]

  # centers <- centers[centers$chr %in% names(coverage),]

  result <- mclapply(names(coverage), function(x) {
    my.cov <- coverage[[x]]
    my.centers <- centers[seqnames(centers) == x, ]
    mw.views <- IRanges::Views(my.cov, start = start(my.centers) -
                                 ceiling(window.size/2), width = window.size + 1)
    flt <- start(mw.views) > 0 & end(mw.views) < length(my.cov)
    mw.views <- mw.views[flt, ]
    my.centers <- my.centers[flt]
    if (length(mw.views) > 0) {
      mat <- as.matrix(mw.views)
      colnames(mat) <- seq(from = (0 - ceiling(window.size/2)),
                           to = 0 + ceiling(window.size/2))
      rownames(mat) <- names(centers)
      return(mat)
    }
    else {
      return(NULL)
    }
  }, mc.cores = n.cores)

  mat <- Reduce(rbind, result)
  centers <- centers[names(centers) %in% rownames(mat),]
  o <- match(names(centers), rownames(mat))
  mat <- mat[o, ]
  if (sum(strand(centers) == "-") > 0) {
    mat[as.vector(strand(centers) == "-"), ] <- t(apply(mat[as.vector(strand(centers) == "-"), ], 1, rev))
  }
  colnames(mat) <- seq(-ceiling(window.size/2), ceiling(window.size/2))
  mat
}


#' Get windowed strand-oriented coverages
#'
#' Out-of-bound windows will be removed!
#'
#' @param windows a named GRanges object (\code{\link[GenomicRanges]{Granges}}
#' @param coverage a coverage object (\code{\link[IRanges]{RleList}} as returned by \code{\link[IRanges]{coverage}})
#'
#' @return a matrix
#'
#' @export
coverageWindowsStranded <- function(windows,  coverage) {

  ## only coverages with enters in centers
  coverage <- coverage[names(coverage) %in% seqlevels(windows)]

  # windows <- windows[seqnames(windows) %in% names(coverage)]

  #  result <- parSapply(cl, names(cov), function(x) {
  result <- lapply(names(coverage), function(x) {
    my.cov <- coverage[[x]]
    my.windows <- windows[seqnames(windows)==x,]
    mw.views <- IRanges::Views(my.cov, start=start(my.windows), end(my.windows))
    ## remove out-of bounds views
    flt <- start(mw.views)>0 & end(mw.views) < length(my.cov)
    mw.views <- mw.views[flt,]
    my.windows <- my.windows[flt,]
    if (length(mw.views) > 0) {
      mat <- as.matrix(mw.views)
      rownames(mat) <- names(my.windows)
      return(mat)
    } else {
      return(NULL)
    }
  })
  #  stopCluster(cl)
  mat <- Reduce(rbind, result)
  windows <- windows[names(windows) %in% rownames(mat),]
  match(names(windows), rownames(mat)) -> o
  mat <- mat[o,]
  if (sum(strand(windows) == "-") > 0) {
    mat[as.vector(strand(windows) == "-"), ] <- t(apply(mat[as.vector(strand(windows) == "-"), ], 1, rev))
  }
  mat
}

#' Get windowed strand-oriented coverages, parallelized
#'
#' Out-of-bound windows will be removed!
#'
#' @param windows a data frame of windows with columns 'chr', 'start', 'end', 'strand'
#' @param coverage a coverage object (\code{\link[IRanges]{RleList}} as returned by \code{\link[IRanges]{coverage}})
#' @param n.cores number of processor cores to use. If not defined all detected cores will be used)
#'
#' @return a matrix
#'
#' @export
coverageWindowsStrandedParallel <- function(windows,  coverage, n.cores=2) {

  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("library 'parallel' is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (is.null(n.cores)) n.cores <- detectCores()


  ## only coverages with enters in centers
  coverage <- coverage[names(coverage) %in% seqlevels(windows)]

  # windows <- windows[seqnames(windows) %in% names(coverage)]

  result <- mclapply(names(coverage), function(x) {
    my.cov <- coverage[[x]]
    my.windows <- windows[seqnames(windows)==x,]
    mw.views <- IRanges::Views(my.cov, start=start(my.windows), end(my.windows))
    ## remove out-of bounds views
    flt <- start(mw.views)>0 & end(mw.views) < length(my.cov)
    mw.views <- mw.views[flt,]
    my.windows <- my.windows[flt,]
    if (length(mw.views) > 0) {
      mat <- as.matrix(mw.views)
      rownames(mat) <- names(my.windows)
      return(mat)
    } else {
      return(NULL)
    }
  }, mc.cores = n.cores)
  mat <- Reduce(rbind, result)
  windows <- windows[names(windows) %in% rownames(mat),]
  match(names(windows), rownames(mat)) -> o
  mat <- mat[o,]
  if (sum(strand(windows) == "-") > 0) {
    mat[as.vector(strand(windows) == "-"), ] <- t(apply(mat[as.vector(strand(windows) == "-"), ], 1, rev))
  }
  mat
}

