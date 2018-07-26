#' Get ocampo paramters per instance
#'
#' @param coverage a coverage object (\code{\link[IRanges]{RleList}} as returned by \code{\link[IRanges]{coverage}})
#' @param references reference regions data frame with columns (chr, start, end, strand). rownames will serve as id for output table
#' @param beforeRef bps before reference to take for shift estimation
#' @param afterRef bps after reference to take for shift estimation
#' @param smoothingWindow size of the smoothing window.
#' @param spacing.low minimum spacing allowed
#' @param spacing.high maximum spacing allowed
#' @param shift.low minimum shift allowed
#' @param shift.high maximum shift allowed
#' @param spacing.high maximum spacing allowed
#'
#' @return a data frame with columns r (correlation), space, and shift for each reference provided
#'
#' @export

ocampo <-
  function(coverage = coverage,
           references = references,
           beforeRef = 200,
           afterRef = 1000,
           smoothingWindow = 100,
           spacing.low = 130,
           spacing.high = 220,
           shift.low = -70,
           shift.high = 60,
           mc.cores = 1) {

    if (!requireNamespace("zoo", quietly = TRUE)) {
      stop("library 'zoo' is needed for this function to work. Please install it.",
           call. = FALSE)
    }
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("library 'parallel' is needed for this function to work. Please install it.",
           call. = FALSE)
    }

    gaussmf <- function (x, sigma, mean) {
      height <- 1
      mfVals = (height * exp(-(((x - mean)^2)/(2 * (sigma^2)))))
    }

    ## generate pattern
    Pattern <- list()
    for (d in spacing.low:spacing.high) {
      Pattern[[d]] <-
        apply(sapply(0:9, function(x) {
          gaussmf(seq(-beforeRef, afterRef), 40, x * d)
        }), 1, sum)
    }

    windows <-
      data.frame(
        chr = references$chr,
        ### the -1+ is a hack to get same vector length after rollmean
        start = ifelse(
          references$strand == "+",
          references$start - (beforeRef + (smoothingWindow / 2)),
          references$end - (afterRef + (smoothingWindow / 2))
        ),
        end = ifelse(
          references$strand == "+",
          references$start + (afterRef + (-1 + smoothingWindow / 2)),
          references$end + (beforeRef + (-1 + smoothingWindow / 2))
        ),
        strand = references$strand
      )
    rownames(windows) <- rownames(references)

    mat <- coverageWindowsStranded(windows, coverage)

    res <-
      parallel::mclapply(1:nrow(mat) , mc.cores = mc.cores , function(ridx) {
        rx <- mat[ridx, ]
        x <- zoo::rollmean(rx, smoothingWindow)
        bestR <- 0
        spacingV <- NA
        shiftV <- NA

        #	if (sum(rx)>10) {
        if (round(sd(x), 1) != 0) {
          ## no reads in window

          for (d in spacing.low:spacing.high) {
            for (shift in shift.low:0) {
              y <- Pattern[[d]]
              r <- cor(x[1:(length(x) + shift)], y[(1 - shift):length(y)])
              if (r > bestR) {
                bestR <- r
                shiftV <- shift
                spacingV <- d
              }
            }
            for (shift in 1:shift.high) {
              y <- Pattern[[d]]
              r <-
                cor(x[(1 + shift):length(x)], y[1:(length(y) - shift)])
              if (r > bestR) {
                bestR <- r
                shiftV <- shift
                spacingV <- d
              }
            }
          }
        }
        c(r=round(bestR, 2), space=spacingV, shift=shiftV)
      })
    df <- t(as.data.frame(res))
    rownames(df) <- rownames(mat)
    df
  }

#' Get ocampo parameters per instance (Tamas version using ccv)
#'
#' @param coverage a coverage object (\code{\link[IRanges]{RleList}} as returned by \code{\link[IRanges]{coverage}})
#' @param references reference regions data frame with columns (chr, start, end, strand). rownames will serve as id for output table
#' @param beforeRef bps before reference to take for shift estimation
#' @param afterRef bps after reference to take for shift estimation
#' @param smoothingWindow size of the smoothing window.
#' @param spacing.low minimum spacing allowed
#' @param spacing.high maximum spacing allowed
#' @param shift.low minimum shift allowed
#' @param shift.high maximum shift allowed
#' @param spacing.high maximum spacing allowed
#' @param sigma_scaled scaling sigma for array model (recommended). Default FALSE
#'
#' @return a data frame with columns r (correlation), space, and shift for each reference provided
#'
#' @export

ocampo2 <- function (coverage = coverage, references = reference, beforeRef = 200,
                     afterRef = 1000, smoothingWindow = 100, spacing.low = 130,
                     spacing.high = 220, shift.low = -70, shift.high = 70, mc.cores = 2,
                     sigma_scaled = FALSE)
{

  gaussmf <- function(x, sigma, mean) {
    height <- 1
    mfVals = (height * exp(-(((x - mean)^2)/(2 * (sigma^2)))))
  }

  Pattern <- list()


  if(sigma_scaled){
    for (d in spacing.low:spacing.high) {
      Pattern[[d]] <- apply(sapply(0:9, function(x) {
        ### !!! scale sigma !!! ###
        gaussmf(seq(-beforeRef, afterRef), 40*(d/150), x * d)
      }), 1, sum)
    }
  } else {
    for (d in spacing.low:spacing.high) {
      Pattern[[d]] <- apply(sapply(0:9, function(x) {
        ### !!! scale sigma !!! ###
        gaussmf(seq(-beforeRef, afterRef), 40, x * d)
      }), 1, sum)
    }
  }

  windows <- data.frame(chr = references$chr,
                        start = ifelse(references$strand == "+",
                                       references$start - (beforeRef + (smoothingWindow/2)),
                                       references$end - (afterRef + (smoothingWindow/2))),
                        end = ifelse(references$strand == "+",
                                     references$start + (afterRef + (-1 + smoothingWindow/2)),
                                     references$end + (beforeRef + (-1 + smoothingWindow/2))),
                        strand = references$strand)
  rownames(windows) <- rownames(references)
  mat <- coverageWindowsStranded(windows, coverage)
  res <- parallel::mclapply(1:nrow(mat), mc.cores = mc.cores,
                            function(ridx) {
                              rx <- mat[ridx, ]
                              x <- zoo::rollmean(rx, smoothingWindow)
                              bestR <- 0
                              spacingV <- NA
                              shiftV <- NA
                              if (round(sd(x), 1) != 0) {
                                for (d in spacing.low:spacing.high) {

                                  #
                                  # This part changed !
                                  #

                                  y <- Pattern[[d]]

                                  my_ccf <- ccf(x,y, lag.max = abs(shift.low), plot=FALSE)

                                  r <- max(my_ccf$acf)
                                  shift <- my_ccf$lag[my_ccf$acf == r]

                                  if (r > bestR) {
                                    bestR <- r
                                    shiftV <- shift
                                    spacingV <- d


                                  }

                                  #
                                  # This part changed !
                                  #
                                }
                              }
                              c(r = round(bestR, 2), space = spacingV, shift = shiftV)
                            })
  df <- t(as.data.frame(res))
  rownames(df) <- rownames(mat)
  df
}








