% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scAPAtrap_funlib.R
\name{extractBcAndUb}
\alias{extractBcAndUb}
\title{Extract barcode and UMI in reads with umi_tools}
\usage{
extractBcAndUb(umitools.path, pattern, input.seq, whitelist, ...)
}
\arguments{
\item{umitools.path}{The path of the umi_tools.}

\item{pattern}{The pattern used to represent barcode and UMI (e.g If the barcode length is 8 and the UMI length is 4, then the pattern=CCCCCCCCNNNN)}

\item{input.seq}{Path of the fastq file, which must contain two paths with the first one recording barcode and UMI information.}

\item{whitelist}{whitelist barcode.}

\item{...}{Arguments passed to other methods and/or advanced arguments.
Advanced arguments:
\describe{
\item{verbose }{If `TRUE` basic status updates will be printed along the way.}
\item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
\item{notRun }{Default is FALSE. If `TRUE`, the Shell commands inside this function will not run but output the commad line.}
}}
}
\value{
The path of Read2 (<input.seq[2]>.extracted.fastq)
}
\description{
Extract barcode and UMI in reads with umi_tools
}
\examples{
umitools.path <-  '/home/dell/miniconda3/envs/scAPA/bin/umi_tools'
input.seq <- c('./fastq/SRR1610598_1.fastq','./fastq/SRR1610598_2.fastq')
pattern <- 'CCCCCCCCNNNN'
whitelist <- './fastq/whitelist.txt'
\dontrun{
extractBcAndUb(umitools.path, pattern, input.seq, whitelist)
}
extractBcAndUb(umitools.path, pattern, input.seq, whitelist, notRun=TRUE)
}
