
#' Get windowed strand-oriented coverages around center points
#'
#' Out-of-bound windows will be removed!
#'
#' @param centers a data frame of center points with columns 'chr', 'center', 'strand'
#' @param window.size the size of the window surrounding the center position, default 1000
#' @param coverage a coverage object (\code{\link[IRanges]{RleList}} as returned by \code{\link[IRanges]{coverage}})
#'
#' @return a matrix
#'
#' @export
coverageWindowsCenteredStranded <- function(centers, window.size=1000, coverage) {

  centers <- centers[centers$chr %in% names(coverage),]

  result <- lapply(names(coverage), function(x) {
    my.cov <- coverage[[x]]
    my.centers <- centers[centers$chr==x,]
    mw.views <- IRanges::Views(my.cov, start=my.centers$center-ceiling(window.size/2), width=window.size+1)
    ## remove out-of bounds views
    flt <- start(mw.views)>0 & end(mw.views) < length(my.cov)
    mw.views <- mw.views[flt,]
    my.centers <- my.centers[flt,]
    if (length(mw.views) > 0) {
      mat <- as.matrix(mw.views)
      colnames(mat) <- seq(from=(0-ceiling(window.size/2)), to=0+ceiling(window.size/2))
      rownames(mat) <- rownames(my.centers)
      return(mat)
    } else {
      return(NULL)
    }
  })

  mat <- Reduce(rbind, result)
  centers <- centers[rownames(centers) %in% rownames(mat),]
  match(rownames(centers), rownames(mat)) -> o
  mat <- mat[o,]
  mat[centers$strand=="-",] <- t(apply(mat[centers$strand=="-",],1,rev))
  colnames(mat) <- seq(-ceiling(window.size/2),ceiling(window.size/2))
  mat
}

#' Get windowed strand-oriented coverages
#'
#' Out-of-bound windows will be removed!
#'
#' @param windows a data frame of windows with columns 'chr', 'start', 'end', 'strand'
#' @param coverage a coverage object (\code{\link[IRanges]{RleList}} as returned by \code{\link[IRanges]{coverage}})
#'
#' @return a matrix
#'
#' @export
coverageWindowsStranded <- function(windows,  coverage) {

  #  cl <- makeCluster(getOption("cl.cores", 8))
  #  clusterExport(cl, list("centers","coverage","window.size") , envir=environment())

  windows <- windows[windows$chr %in% names(coverage),]

  #  result <- parSapply(cl, names(cov), function(x) {
  result <- lapply(names(coverage), function(x) {
    my.cov <- coverage[[x]]
    my.windows <- windows[windows$chr==x,]
    mw.views <- IRanges::Views(my.cov, start=my.windows$start, my.windows$end)
    ## remove out-of bounds views
    flt <- start(mw.views)>0 & end(mw.views) < length(my.cov)
    mw.views <- mw.views[flt,]
    my.windows <- my.windows[flt,]
    if (length(mw.views) > 0) {
      mat <- as.matrix(mw.views)
      rownames(mat) <- rownames(my.windows)
      return(mat)
    } else {
      return(NULL)
    }
  })
  #  stopCluster(cl)
  mat <- Reduce(rbind, result)
  windows <- windows[rownames(windows) %in% rownames(mat),]
  match(rownames(windows), rownames(mat)) -> o
  mat <- mat[o,]
  mat[windows$strand=="-",] <- t(apply(mat[windows$strand=="-",],1,rev))
  mat
}
