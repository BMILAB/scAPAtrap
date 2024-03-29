% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scAPAtrap_funlib.R
\name{findTails}
\alias{findTails}
\title{Find the positions with polyA tails genome/chromosome-wide}
\usage{
findTails(bamfile, chr = NULL, len = NULL, tailsfile = NULL, ...)
}
\arguments{
\item{bamfile}{A BAM file which needs the index file (.bai).}

\item{chr}{NULL or a character string. If NULL, search whole genome. Otherwise search tails on one chr.}

\item{len}{NULL or the chromosome length of the given chr.}

\item{tailsfile}{if not NULL then output to tailsfile, and return the tailsfile name.}

\item{...}{Arguments passed to other methods and/or advanced arguments.
Advanced arguments:
\describe{
\item{verbose }{If `TRUE` basic status updates will be printed along the way.}
\item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
}}
}
\value{
A data frame recording tail positions with chr/strand/coord/count; or tailsfile name (if tailsfile not NULL)
}
\description{
Find the precise positions with polyA tails genome-wide or on one given chromosome.
             Given chr, this function will traverse each chromosome and find the site through the A-rich at the end.
}
\examples{
\dontrun{
findTails(bamfile, chr=1, len=1000000001)
findTails(bamfile, tailsfile='output/output.tails')
}
}
