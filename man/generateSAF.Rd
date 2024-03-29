% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scAPAtrap_funlib.R
\name{generateSAF}
\alias{generateSAF}
\title{Generate a peak annotation file}
\usage{
generateSAF(forwardPeaks, reversePeaks, outputdir, ...)
}
\arguments{
\item{forwardPeaks}{A peak filename or a peak data.frame of forward strand, which is the output of \code{\link{findPeaksByStrand}} or \code{\link{findPeaks}}, where strand = +.}

\item{reversePeaks}{Save as forwardPeaks but on the reverse strand. Can provide forwardPeaks or reversePeaks, or both.}

\item{outputdir}{Output file directory.}

\item{...}{Arguments passed to other methods and/or advanced arguments.
Advanced arguments:
\describe{
\item{verbose }{If `TRUE` basic status updates will be printed along the way.}
\item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
}}
}
\value{
Full path of the peak list (<outputdir>/peaks.saf). This file contains five columns without header: peakID/chr/start/end/strand.
}
\description{
Generate a peak annotation file (PeakID/chr/start/end/Strand) by combining forward peaks and reverse peaks generated by findPeaks
}
\examples{
\dontrun{
forwardPeaks <-findPeaks(fullcovF,'+', 98, 1000, cutoff = 1)
reversePeaks <-findPeaks(fullcovR,'-',98, 1000, cutoff = 1)
generateSAF(forwardPeaks, reversePeaks,'./data')
}
}
