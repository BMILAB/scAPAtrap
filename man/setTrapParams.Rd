% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scAPAtrap_funlib.R
\name{setTrapParams}
\alias{setTrapParams}
\title{Set the parameters for running scAPAtrap}
\usage{
setTrapParams(trap.params = NULL, print = TRUE, check = TRUE, ...)
}
\arguments{
\item{trap.params}{A list similar to \code{\link{TRAP.PARAMS}}.}

\item{print}{Whether to print the list (default is TRUE).}

\item{check}{Whether to check the validity of values (default is TRUE).}

\item{...}{Arguments passed to other methods and/or advanced arguments.
Advanced arguments:
\describe{
\item{verbose }{If `TRUE` basic status updates will be printed along the way.}
\item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
}}
}
\value{
A list like code{\link{TRAP.PARAMS}}.
}
\description{
setTrapParams sets parameters for running scAPAtrap, by setting default values using \code{\link{TRAP.PARAMS}} and checking the validity.
}
\examples{
## get default
trap.params=setTrapParams()
## change some parameters
trap.params=setTrapParams(list(chrs=c('2L','2R','3L','3R','4','X','Y'), TenX=FALSE))
## will ignore non-valid parameter names
trap.params=setTrapParams(list(xx=FALSE))
## add barcodes
trap.params=setTrapParams(list(barcode=c('AAA','BBB')))
## not print
trap.params=setTrapParams(print=FALSE)
}
