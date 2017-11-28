
#' Cumulative plot with interquartile ranges
#'
#' @param mat matrix or data.frame to cumulate. colum names specify x
#' @param base.col base color
#'
#' @export

cumPlot <- function(mat, base.col="#000000",...) {
  x <- as.integer(colnames(mat))
  y <- apply(mat,2,mean)
  y1 <- apply(mat,2,quantile, probs=0.25)
  y3 <- apply(mat,2,quantile, probs=0.75)
  ylim <- c(0,max(c(y, y1, y3)))
  plot(x, y, type="n", ylim=ylim, xlab="distance to TSS (bp)", ylab="dyad density")
  polygon(c(x,rev(x)),c(y3,rev(y1)),col=paste(base.col,"20", sep=""), border=paste(base.col,"50", sep=""))
  lines(x, y, lwd=2, col=base.col)
}
