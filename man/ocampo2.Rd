% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ocampo.R
\name{ocampo2}
\alias{ocampo2}
\title{Get ocampo parameters per instance (Tamas version using ccv)}
\usage{
ocampo2(coverage = coverage, references = reference, beforeRef = 200,
  afterRef = 1000, smoothingWindow = 100, spacing.low = 130,
  spacing.high = 220, shift.low = -70, shift.high = 70, mc.cores = 2,
  sigma_scaled = FALSE)
}
\arguments{
\item{coverage}{a coverage object (\code{\link[IRanges]{RleList}} as returned by \code{\link[IRanges]{coverage}})}

\item{references}{reference regions data frame with columns (chr, start, end, strand). rownames will serve as id for output table}

\item{beforeRef}{bps before reference to take for shift estimation}

\item{afterRef}{bps after reference to take for shift estimation}

\item{smoothingWindow}{size of the smoothing window.}

\item{spacing.low}{minimum spacing allowed}

\item{spacing.high}{maximum spacing allowed}

\item{shift.low}{minimum shift allowed}

\item{shift.high}{maximum shift allowed}

\item{sigma_scaled}{scaling sigma for array model (recommended). Default FALSE}

\item{spacing.high}{maximum spacing allowed}
}
\value{
a data frame with columns r (correlation), space, and shift for each reference provided
}
\description{
Get ocampo parameters per instance (Tamas version using ccv)
}
