#' @import GenomicRanges
#' @importFrom dplyr group_by summarise mutate
#' @importFrom Matrix sparseMatrix
#' @importFrom reshape2 dcast
#' @importFrom IRanges IRanges
#' @importFrom utils write.table
#' @importFrom stats quantile
#' @importFrom magrittr %>%
NULL

## update LOG (2023-09-10):
## Fix a bug in findPeaks: previously, when no peaks wider than maxwidth, will not output small peaks.
## Add function findPeaksByStrand, which combines previous loadBpCoverages and findPeaks.
## Add wrapper function scAPAtrap for one step running.
## Add initScAPAtrap and demo script for step by step run of scAPAtrap.
## Add neat LOG and message.
## Add function findTails (genome-wide searching) and findTailsByPeaks (peaks-wide searching) to find tails independently.
## Add function reducePeaks to remove lowly expressed peaks, which can be called before searching tails or output.
## Add function convertAPAtrapData to convert scAPAtrapData (a list) to other objects incuding PACdataset, SeurateObject, SingleCellExperiment.
## Update function generatescExpMa to allow tail=NULL for not searching tails of peaks; output a list containing sparseMatrix for bigger data.
## More flexible parameter settings (TRAP.PARAMS): allowing flexible chrs, barcode, tails.
## Add global variable TRAPFILES to record output files in scAPAtrap.

## ??
## Different runs of scAPAtrap generate a very tiny difference, random seed?

##---- build --------
## rm(list=ls())
## devtools::check()
## devtools::load_all()
## roxygen2::roxygenise()
## >>> Writing LOG.Rd
## build&reload

#library(conflicted)
# library(dplyr)
# library(reshape2)
# library(Matrix)
# library(GenomicRanges)

########## Utils ####################
.getDotArg<-function (name, value, ...) {
  args <- list(...)
  if (!name %in% names(args)) {
    return(value)
  }
  else {
    return(args[[name]])
  }
}

## global variable to record output files in scAPAtrap, a two-column data.frame (file, what)
## files' what marked with * are temporary files that can be savely deleted.
myenv=new.env()
myenv$TRAPFILES <- data.frame(matrix(nrow=0, ncol=2, dimnames=(list(c(), c('file','what')))))
#utils::globalVariables(c("TRAPFILES", "TRAP.PARAMS")) #not works

#' TRAPFILES global variable
#'
#' TRAPFILES is a global variable after running scAPAtrap, containing two columns (file, what).
#' This variable is automatically updated when running certain steps in scAPAtrap.
#' This function get TRAPFILES or clear TRAPFILES, but cannot add/modify items (rows) in TRAPFILES.
#' @param clear TRUE to clear TRAPFILES variable
#' @return A data frame of files containing two columns (file, what).
#' @examples
#' TRAPFILES()
#' countPeaks(umitools.path='umitools.path', input='xx.bam', outputdir='xxdir', TenX=TRUE, notRun=TRUE)
#' TRAPFILES()
#' TRAPFILES(clear=TRUE)
#' @export
TRAPFILES<-function(clear=FALSE) {
  if (clear) {
    myenv$TRAPFILES=myenv$TRAPFILES[character(0), ]
  }
  return(myenv$TRAPFILES)
}

## add file name and description in each function, those files marked with * can be savely deleted
.addTrapFiles<-function(files, whats) {
  if (length(files)>length(whats)) whats=rep('unknown', 1:length(files))
  for (i in 1:length(files)) {
    if (!(files[i] %in% myenv$TRAPFILES$file)) myenv$TRAPFILES <- rbind(myenv$TRAPFILES, c(files[i], whats[i]))
    colnames(myenv$TRAPFILES) <- c('file','what')
  }
}

## output a string to logf or screen
## It is easily to output a string with >>> and time to the screen and/or logfile
# .msg('command ...', verbose=FALSE) ## nothing
# .msg('command ...') --> >>> 2023-09-06 16:48:15 generatescExpMa: command...
# .msg('command ...', verbose=TRUE, logf='xx.log') --> >>> 2023-09-06 16:48:15 generatescExpMa: command...
# .msg('##', pre='', showTime=FALSE) --> ##
.msg<-function(string='', pre='>>>', showTime=TRUE, newLine=TRUE, preLine=FALSE, ...) {
  verbose <- .getDotArg("verbose", TRUE, ...)
  logf <- .getDotArg("logf", NULL, ...)

  t=ifelse(showTime, paste0(as.character(Sys.time()),' '), '')
  nl=ifelse(newLine, '\n', '')
  pl=ifelse(preLine, '\n', '')
  if (pre!='') pre=paste0(pre, ' ')
  s=sprintf("%s%s%s%s%s", pl, pre, t, string, nl)
  if (verbose) cat(s)
  if (is.character(logf)) {
    cat(s, file=logf, append=TRUE)
  }
}

# show msg within functions, by adding ... + function name + string.
# .msgin(fn='function', string='command***', verbose=TRUE, logf='xx.log') --> ... 2023-09-06 17:13:31 function: command***
# .msgin(fn='', string='command***', verbose=TRUE, logf=NULL) --> ... 2023-09-06 17:21:32 command***
# .msgin(fn='FUNC', string='command***', verbose=TRUE, logf=NULL) --> ... 2023-09-06 17:21:52 FUNC: command***
.msgin<-function(fn='', string='', ...) {
  if (fn=='') s=string else s=paste0(fn, ': ', string)
  .msg(string=s, pre='...', showTime=TRUE, newLine=TRUE, ...)
}

## run command line and output message on screen or file
## verbose=TRUE to show only the command line; FALSE now show
## logf=file to sink command and command output to file but NOT show on the screen; logf=NULL to show command output on the screen
## however, system() not allows capture screen output of inner command, the commond message will be shown on screen no matter verbose=FALSE
## .runCommand('', command="tree") # ... 2023-09-06 18:01:22 tree +++ tree results on screen +++ ... 2023-09-06 18:01:22 command done.
## .runCommand('', command="tree", logf='xx.log')  # ... 2023-09-06 18:01:22 tree +++ tree results in xx.log +++ ... 2023-09-06 18:01:22 command done.
## .runCommand('', command="tree", verbose=FALSE, logf='xx.log') # nothing on screen, but tree results in xx.log
## .runCommand('', command="tree", verbose=FALSE, logf=FALSE) # nothing on screen or file
.runCommand<-function(fn='', command,  ofile=NULL, ocheck=FALSE, ...) {

  verbose <- .getDotArg("verbose", TRUE, ...)
  notRun <- .getDotArg("notRun", FALSE, ...)
  wait <- .getDotArg("wait", TRUE, ...)
  logf <- .getDotArg("logf", TRUE, ...)

  if (notRun & !verbose &is.null(logf)) verbose=TRUE

  .msgin(fn=fn, string=command, ...)

  ## not run, only output command line string
  if (notRun) return(invisible(NULL))

  ## run command
  #s=system2(command = command, wait = wait, stdout=TRUE, stderr=TRUE) # return str vectors
  s=system(command = command, wait = wait, intern=FALSE, ignore.stdout=FALSE,  ignore.stderr=FALSE)

  if (s!=0) {
    warning("It seems the command failed (exit status!=0), please see screen message or log file for details!\n")
  }
  if (ocheck) {
    if (!is.null(ofile)) {
      if (!file.exists(ofile)) stop(ofile, "not exists! Probably the command failed, please see screen message or log file for detailes!")
    }
  }

  # if (length(s)>0) {
  #   if (is.character(logf)) { #output to file, now show on screen
  #     cat(s, file = logf, append=TRUE, sep='\n')
  #   } else {# on screen
  #     cat(s, sep='\n')
  #   }
  # }
  .msgin(fn=fn, 'command done.', ...)
}

#-----------------------------------Preprocess----------------------------------------------
#' Filter a BAM file to preserve unique mappings with samtools
#'
#' findUniqueMap filters a BAM file to preserve unique mappings with samtools. This function will generate .bam (<input>.uniq.bam or UniqSorted.bam if sort=TRUE) and .bai (if index=TRUE).
#' Suffix .Uniq or .UniqSorted will be added to the input BAM file. This function is the same as call samtools in Shell.
#'
#' @param samtools.path The path of the samtools.
#' @param input BAM file name.
#' @param thread Number of CPU threads, default is 12.
#' @param sort Logical value, TRUE to sort the BAM file.
#' @param index Logical value, TRUE to build the index.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{If `TRUE` basic status updates will be printed along the way.}
#' \item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
#' \item{notRun }{Default is FALSE. If `TRUE`, the Shell commands inside this function will not run but output the commad line.}
#' }
#' @examples
#' samtools.path <- '/home/aa/miniconda2/envs/nar_env/bin/samtools'
#' input <- './data/demo.bam'
#' \dontrun{
#' findUniqueMap(samtools.path, input, 24)
#' }
#' ## only output command lines
#' findUniqueMap(samtools.path, input, 24, notRun=TRUE)
#' @return Output path of the new BAM file name with only unique mapped reads (<input>.uniq.bam or UniqSorted.bam if sort=TRUE).
#' @export
findUniqueMap <- function(samtools.path, input, thread=12, sort = TRUE, index = TRUE, ...) {

  fn='findUniqueMap'

  .msg(sprintf("%s: start.", fn), preLine=FALSE, ...)

  output <- paste0(gsub(".bam$", "", input), ".Uniq.bam")
  command <- paste0(samtools.path, " view -@ ", thread, " -h -F 256 -bS ", input, " > ", output)
  .runCommand(fn, command, ofile=output, ocheck=TRUE, ...)
  #.addTrapFiles(output, 'findUniqueMap: samtools view unique')

  res <- output
  if (sort) {
      output.sorted <- paste0(gsub(".Uniq.bam$", "", output), ".UniqSorted.bam")
      command <- paste0(samtools.path, " sort -@ ", thread, " -o ", output.sorted, " ", output)
      .runCommand(fn, command, ofile=output.sorted, ocheck=TRUE, ...)
      unlink(res) # delete uniq.bam if sorted
      res <- output.sorted
      .addTrapFiles(res, '*findUniqueMap: samtools sort')
  }
  if (index) {
      command <- paste0(samtools.path, " index -@ ", thread, " ", res)
      .runCommand(fn, command, ...)
      .addTrapFiles(paste0(res,'.bai'), '*findUniqueMap: samtools index')
  }

  .msg(sprintf("%s: finish.", fn), ...)

  return(res)
}


#' Deduplication with umitools
#'
#' @param umitools.path The path of the umi_tools.
#' @param input A BAM file which needs the index file.
#' @param TenX Logical value, TRUE for 10X data.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{If `TRUE` basic status updates will be printed along the way.}
#' \item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
#' \item{notRun }{Default is FALSE. If `TRUE`, the Shell commands inside this function will not run but output the commad line.}
#' }
#' @examples
#' umitools.path <- '/home/aa/miniconda2/envs/umi_tools/bin/umi_tools'
#' input <- './data/demo.Uniq.sorted.bam'
#' \dontrun{
#' dedupByPos(umitools.path, input)
#' }
#' dedupByPos(umitools.path, input, notRun=TRUE)
#' @return Output path of the deduplicated BAM file (<input>.dedup.bam).
#' @export
dedupByPos <- function(umitools.path, input, TenX = TRUE, ...) {

  fn='dedupByPos'
  .msg(sprintf("%s: start.", fn), preLine=FALSE, ...)

  output <- paste0(gsub(".bam$", "", input), ".dedup.bam")
  if (TenX) {
    command <- paste0(umitools.path, " dedup -I ", input, " -S ", output, " --method=unique --extract-umi-method=tag --umi-tag=UB --cell-tag=CB")
  } else {
    command <- paste0(umitools.path, " dedup -I ", input, " -S ", output, " --method=unique")
  }

  .runCommand(fn, command, ofile=output, ocheck=TRUE, ...)

  .msg(sprintf("%s: finish.", fn), ...)

  .addTrapFiles(output, '*dedupByPos: umitools dedup')
  return(output)
}


#' Separate BAM files by strand
#'
#'@param samtools.path The path of the samtools.
#'@param input The BAM file name.
#'@param thread Number of CPU threads, default is 12.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{If `TRUE` basic status updates will be printed along the way.}
#' \item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
#' \item{notRun }{Default is FALSE. If `TRUE`, the Shell commands inside this function will not run but output the commad line.}
#' }
#'@examples
#' samtools.path <- '/home/aa/miniconda2/envs/nar_env/bin/samtools'
#' input <- './data/demo.Uniq.sorted.dedup.bam'
#' \dontrun{
#' separateBamBystrand(samtools.path,input,24)
#' }
#' separateBamBystrand(samtools.path,input,24, notRun=TRUE)
#'@return String vector storing the path of the forward BAM file (<input>.forward.bam) and the reverse BAM file (<input>.reverse.bam)
#'@export
separateBamBystrand <- function(samtools.path, input, thread=12, ...) {

  fn='separateBamBystrand'

  .msg(sprintf("%s: start.", fn), preLine=FALSE, ...)

  outputF <- paste0(gsub(".bam", "", input), ".forward.bam")
  command <- paste0(samtools.path, " view -@ ", thread, " -h -F 0x10 -bS ", input, " > ", outputF)
  .runCommand(fn, command, ofile=outputF, ocheck=TRUE, ...)
  .addTrapFiles(outputF, '*separateBamBystrand: samtools forward BAM')

  outputR <- paste0(gsub(".bam", "", input), ".reverse.bam")
  command <- paste0(samtools.path, " view -@ ", thread, " -h -f 0x10 -bS ", input, " > ", outputR)
  .runCommand(fn, command, ofile=outputR, ocheck=TRUE,  ...)
  .addTrapFiles(outputR, '*separateBamBystrand: samtools reverse BAM')

  command <- paste0(samtools.path, " index -@ ", thread, " ", outputF)
  .runCommand(fn, command, ...)
  .addTrapFiles(paste0(outputF,'.bai'), '*separateBamBystrand: samtools forward index')

  command <- paste0(samtools.path, " index -@ ", thread, " ", outputR)
  .runCommand(fn, command, ...)
  .addTrapFiles(paste0(outputR,'.bai'), '*separateBamBystrand: samtools reverse index')

  .msg(sprintf("%s: finish.", fn), ...)

  return(c(outputF, outputR))
}


#' Extract barcode and UMI in reads with umi_tools
#'
#' @param umitools.path The path of the umi_tools.
#' @param pattern The pattern used to represent barcode and UMI (e.g If the barcode length is 8 and the UMI length is 4, then the pattern=CCCCCCCCNNNN)
#' @param input.seq Path of the fastq file, which must contain two paths with the first one recording barcode and UMI information.
#' @param whitelist whitelist barcode.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{If `TRUE` basic status updates will be printed along the way.}
#' \item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
#' \item{notRun }{Default is FALSE. If `TRUE`, the Shell commands inside this function will not run but output the commad line.}
#' }
#' @examples
#' umitools.path <-  '/home/dell/miniconda3/envs/scAPA/bin/umi_tools'
#' input.seq <- c('./fastq/SRR1610598_1.fastq','./fastq/SRR1610598_2.fastq')
#' pattern <- 'CCCCCCCCNNNN'
#' whitelist <- './fastq/whitelist.txt'
#' \dontrun{
#' extractBcAndUb(umitools.path, pattern, input.seq, whitelist)
#' }
#' extractBcAndUb(umitools.path, pattern, input.seq, whitelist, notRun=TRUE)
#'@return The path of Read2 (<input.seq[2]>.extracted.fastq)
#' @export
extractBcAndUb <- function(umitools.path, pattern, input.seq, whitelist, ...) {

  fn='extractBcAndUb'
  .msg(sprintf("%s: start.", fn), preLine=FALSE, ...)

  output1 <- paste0(gsub(".fastq$", "", input.seq[1]), ".extracted.fastq")
  output2 <- paste0(gsub(".fastq$", "", input.seq[2]), ".extracted.fastq")
  command <- paste0(umitools.path, " extract --bc-pattern=", pattern, " --stdin ", input.seq[1], " --stdout ", output1,
                    " --read2-in ", input.seq[2], " --read2-out=", output2, " --filter-cell-barcode ", " --whitelist=", whitelist)

  .runCommand(fn, command, ofile=output2, ocheck=TRUE, ...)
  .addTrapFiles(output1, '*extractBcAndUb: umitools bc R1')
  .addTrapFiles(output2, '*extractBcAndUb: umitools bc R2')

  .msg(sprintf("%s: finish.", fn), ...)

  return(output2)
}


#' Build the reference genome index with STAR
#'
#' @param star.path The path of the STAR.
#' @param genome.fasta The reference genome sequence file in fasta format.
#' @param indexdir Directory of the reference genome index.
#' @param thread Number of CPU threads, default is 12.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{If `TRUE` basic status updates will be printed along the way.}
#' \item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
#' \item{notRun }{Default is FALSE. If `TRUE`, the Shell commands inside this function will not run but output the commad line.}
#' }
#' @return Directory of reference genome index
#' @examples
#' star.path='/home/dell/miniconda3/envs/scAPA/bin/STAR'
#' genome.fasta <- '/home/aa/genome_data/mm10/Mus_musculus.GRCm38.dna.primary_assembly.fa'
#' indexdir <- './mm10_index/'
#' \dontrun{
#' generateRefIndex(star.path, genome.fasta, indexdir, 24)
#' }
#' generateRefIndex(star.path, genome.fasta, indexdir, 24, notRun=TRUE)
#' @export
generateRefIndex <- function(star.path, genome.fasta, indexdir, thread=12, ...) {


  fn='generateRefIndex'
  .msg(sprintf("%s: start.", fn), preLine=FALSE, ...)

  command <- paste0(star.path, " --runMode genomeGenerate --genomeDir ", indexdir, " --genomeFastaFiles ", genome.fasta,
                    " --runThreadN ", thread)
  .runCommand(fn, command, ofile=indexdir, ocheck=TRUE, ...)
  .addTrapFiles(indexdir, '*generateRefIndex: STAR index genome')

  .msg(sprintf("%s: finish.", fn), ...)
  return(indexdir)
}


#' Align the sequence to the reference genome with STAR
#'
#' @param star.path The path of the STAR.
#' @param indexdir Directory of reference genome index.
#' @param input.seq Input fastq file.
#' @param prefix The prefix of the output file name.
#' @param thread Number of CPU threads, default is 12.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{If `TRUE` basic status updates will be printed along the way.}
#' \item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
#' \item{notRun }{Default is FALSE. If `TRUE`, the Shell commands inside this function will not run but output the commad line.}
#' }
#' @examples
#' star.path='/home/dell/miniconda3/envs/scAPA/bin/STAR'
#' indexdir <- "./"
#' input.seq <- 'SRR1610598_2.extracted.fastq'
#' \dontrun{
#' generateAlignBam(star.path, indexdir, input.seq, prefix='./BAM/rep1',12)
#' }
#' generateAlignBam(star.path, indexdir, input.seq, prefix='./BAM/rep1', notRun=TRUE)
#' @return NULL
#' @export
generateAlignBam <- function(star.path, indexdir, input.seq, prefix, thread=12, ...) {

  fn='generateAlignBam'
  .msg(sprintf("%s: start.", fn), preLine=FALSE, ...)

  command <- paste0(star.path, " --runThreadN ", thread, " --genomeDir ", indexdir, " --readFilesIn ", input.seq, " --outFilterMultimapNmax 1 ",
                    " --outSAMtype BAM SortedByCoordinate ", " --outFileNamePrefix ", prefix)
  .runCommand(fn, command, ...)

  .msg(sprintf("%s: finish.", fn), ...)
  ##bug: output?
}

#-----------------------------------findTails----------------------------------------------
.readChrInfo <- function(file) {
  chr_length <- utils::read.delim(file = file)
  colnames(chr_length) <- c("flag", "chr", "coord")
  chr_length$chr <- gsub("SN:", "", chr_length$chr)
  chr_length$coord <- as.numeric(gsub("LN:", "", chr_length$coord))
  chr_length <- subset(chr_length, chr != "MT")
  chr.g <- with(chr_length, GRanges(seqnames = chr, ranges = IRanges::IRanges(start = 1, end = coord)))

  return(chr.g)
}

## search BAM to find tails, returning chr/strand/coord/count ranges
.findTails<-function(bamfile, bamParam) {
  gal1 <- GenomicAlignments::readGAlignments(bamfile, use.names = TRUE, param = bamParam)
  s_1 <- (grepl("[0-9]*M[1-9]{2,}S", gal1@cigar) & as.vector(gal1@strand) == "+")
  s_2 <- (grepl("[1-9]{2,}S[0-9]*M", gal1@cigar) & as.vector(gal1@strand) == "-")

  bam1 <- gal1[s_1]
  bam2 <- gal1[s_2]

  bam1 <- bam1[grepl("(A{3,}[^A]{0,2})*A{6,}([^A]{0,2}A{3,})*.{0,2}?$", bam1@elementMetadata@listData$seq)]
  bam2 <- bam2[grepl("^.{0,2}?(T{3,}[^T]{0,2})*T{6,}([^T]{0,2}T{3,})*", bam2@elementMetadata@listData$seq)]

  final_bam1 <- data.frame(chr = as.vector(seqnames(bam1)), strand = as.vector(strand(bam1)), coord = end(bam1))
  final_bam2 <- data.frame(chr = as.vector(seqnames(bam2)), strand = as.vector(strand(bam2)), coord = start(bam2))

  bam <- rbind(final_bam1, final_bam2)
  bam <- dplyr::group_by(bam, chr, strand, coord) %>% dplyr::summarise(count = dplyr::n())
  ## warning.. override using .group argument..
  return(bam)
}

.getBAMchrs<-function(bamfile) {
  chrs <- as.character(seqnames(seqinfo(Rsamtools::BamFile(bamfile))))
  return(chrs)
}

.indexBam<-function(bamfile, ...) {
  baifile=paste0(bamfile, '.bai')
  if(!file.exists(baifile)) {
    .msgin(fn='indexBam', bamfile, ...)
    Rsamtools::indexBam(bamfile)
  }
  return(baifile)
}

#' Find the positions with polyA tails genome/chromosome-wide
#'
#' @description Find the precise positions with polyA tails genome-wide or on one given chromosome.
#'              Given chr, this function will traverse each chromosome and find the site through the A-rich at the end.
#' @param bamfile A BAM file which needs the index file (.bai).
#' @param chr NULL or a character string. If NULL, search whole genome. Otherwise search tails on one chr.
#' @param len NULL or the chromosome length of the given chr.
#' @param tailsfile if not NULL then output to tailsfile, and return the tailsfile name.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{If `TRUE` basic status updates will be printed along the way.}
#' \item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
#' }
#' @return A data frame recording tail positions with chr/strand/coord/count; or tailsfile name (if tailsfile not NULL)
#' @examples
#' \dontrun{
#' findTails(bamfile, chr=1, len=1000000001)
#' findTails(bamfile, tailsfile='output/output.tails')
#' }
#' @export
findTails <- function(bamfile, chr=NULL, len=NULL, tailsfile=NULL, ...) {

  fn='findTails'
  .msg(sprintf("%s: start.", fn), preLine=FALSE, ...)

  what <- c("rname", "strand", "pos", "cigar", "seq")
  if (is.null(chr)) {
    param <- Rsamtools::ScanBamParam(what = what)
  } else {

    if (is.null(len)) stop("findTails: chr is provided, then len (the chromosome length) should be also provided!")
    # check the chrnames in bam to be consistent with chr
    # chrs=.getBAMchrs(bamfile)
    # if (!(chr %in% chrs)) stop(sprintf("findTails: chr=%s not in chrs of the bamfile (%s...)!", chr, .shortStr(chrs)))
    .indexBam(bamfile, ...)

    which <- GenomicRanges::GRanges(seqnames = chr, IRanges::IRanges(start=1, end=len))
    param <- Rsamtools::ScanBamParam(what = what, which = which)
  }
  #chr <- colnames(bam)[1]

  res=.findTails(bamfile, param)
  if (!is.null(tailsfile)) {
    .msgin(fn, sprintf("save %d tails-coords to %s", nrow(res), tailsfile), ...)
    write.table(res, file=tailsfile, col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)
    .addTrapFiles(tailsfile, 'findTails: tails-coords')
    .msg(sprintf("%s: finish.", fn), ...)
    return(tailsfile)
  } else {
    .msgin(fn, sprintf("find %d tails-coords", nrow(res)), ...)
    .msg(sprintf("%s: finish.", fn), ...)
    return(res)
  }

}


#' Find tail positions around/within peak ranges
#'
#' @description Find the precise positions with polyA tails within d nt of peak ranges, given a peaksfile recording peak ranges.
#'              This function narrows the search range of tails to peak ranges to speed up tail searching, which is useful for large BAM file.
#' @param bamfile A BAM file which needs the index file (.bai).
#' @param peaksfile The path of the peaks.saf file generated by \code{\link{generateSAF}}.
#' @param d Add a 5' and 3' margin of d nt to each peak to search tails. Default is 200.
#' @param tailsfile if not NULL then output to tailsfile, and return the tailsfile name.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{If `TRUE` basic status updates will be printed along the way.}
#' \item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
#' }
#' @return A data frame recording tail positions with chr/strand/coord/count; or tailsfile name (if tailsfile not NULL).
#' @examples
#' \dontrun{
#' findTailsByPeaks(bamfile, peaksfile, d=500)
#' findTailsByPeaks(bamfile, peaksfile, d=500, tailsfile='output.tails')
#' }
#' @export
findTailsByPeaks <- function(bamfile, peaksfile, d=200, tailsfile=NULL, ...) {

  fn='findTailsByPeaks'
  .msg(sprintf("%s: start.", fn), preLine=FALSE, ...)

  peaks <- utils::read.delim(peaksfile, header = FALSE)
  names(peaks) <- c("peakID", "chr", "start", "end", "strand")
  gr <- with(peaks, GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start = start, end = end), strand = strand))
  if (d<=0) d=0
  if (d>0) {
    gr=GenomicRanges::resize(gr, width=width(gr)+d, fix="end", use.names=TRUE, ignore.strand=FALSE)
    gr=GenomicRanges::resize(gr, width=width(gr)+d, fix="start", use.names=TRUE, ignore.strand=FALSE)
    gr=GenomicRanges::reduce(gr)
  }

  .indexBam(bamfile, ...)

  # validate chr consistency
  cb=.getBAMchrs(bamfile)
  cp=unique(as.character(peaks$chr))
  if (!all(cp %in% cb)) {
    if (!any(cp %in% cb))
      stop(sprintf("findTailsByPeaks: peak chrs (%s...) not in bam chrs (%s...)!", .shortStr(cp), .shortStr(cb)))
  }

  what <- c("rname", "strand", "pos", "cigar", "seq")
  param <- Rsamtools::ScanBamParam(what = what, which = gr)
  .msgin(fn, sprintf("findTailsByPeaks on %d reduced peak ranges (d=%dnt)",length(gr), d), ...)

  res=.findTails(bamfile, param)

  if (!is.null(tailsfile)) {
    .msgin(fn, sprintf("save %d tails-coords to %s", nrow(res), tailsfile), ...)
    write.table(res, file=tailsfile, col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)
    .addTrapFiles(tailsfile, 'findTails: tails-coords')
    .msg(sprintf("%s: finish.", fn), ...)
    return(tailsfile)
  } else {
    .msgin(fn, sprintf("find %d tails-coords", nrow(res)), ...)
    .msg(sprintf("%s: finish.", fn), ...)
    return(res)
  }
}


#-----------------------------------findPeaks----------------------------------------------
#' Load a BAM file and calculate coverages of each position with derfinder
#'
#' @param files BAM files which need the corresponding index files (.bai).
#' @param chrs Chromosome information.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{If `TRUE` basic status updates will be printed along the way.}
#' \item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
#' }
#' @return A list with one element per chromosome. Each element is a DataFrame with the coverage information.
#' @examples
#' \dontrun{
#' chrs <- c(as.character(1:19),'X','Y')
#' fullcovF <- loadBpCoverages('./dedup_h2.forward.sorted.bam', chrs)
#' }
#' @export
loadBpCoverages <- function(files, chrs, ...) {

  fn='loadBpCoverages'
  .msg(sprintf("%s: start.", fn), preLine=FALSE, ...)

  fullCov <- derfinder::fullCoverage(files = files, chrs = chrs, ...)

  ## no SQ lines in @header
  # # validate chr consistency
  # if (is.null(chrs))
  #   chrs=.getBAMchrs(files[1])
  # else if (length(chrs)==0) {
  #   chrs=.getBAMchrs(files[1])
  # } else {
  #   cb=.getBAMchrs(files[1])
  #   if (!all(chrs %in% cb)) {
  #     if (!any(chrs %in% cb))
  #       stop(sprintf("loadBpCoverages: chrs (%s...) not in bam chrs (%s...)!", substr(toString(chrs), 1, 20), substr(toString(cb),1, 20)))
  #   }
  # }
  .msg(sprintf("%s: finish.", fn), ...)

  return(fullCov)
}


#' Perform peak calling on a BAM file
#'
#' Perform peak calling on a BAM file with single strand (+ or -) using the output of loadBpCoverages. It is recommended to use findPeaksByStrand instead.
#'
#' @param fullCov The output of loadBpCoverages.
#' @param strand strand.
#' @param L Read length, used to filter peaks wider than L.
#' @param maxwidth Maximum peak width, used to filter peaks narrower than maxwidth. Wider peaks will be splitted to smaller ones.
#' @param cutoff The base level coverage cutoff to filter valid peaks with enough coverages, default is 10.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{If `TRUE` basic status updates will be printed along the way.}
#' \item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
#' }
#' @return A data frame with five columns: chr/start/end/strand/value, with each row being a peak. The value column is the average coverage of peak.
#' @examples
#' \dontrun{
#' forwardPeaks <-findPeaks(fullcovF, '+', 98, 1000)
#' reversePeaks <-findPeaks(fullcovR, '-', 98, 1000)
#' }
#' @export
findPeaks <- function(fullCov, strand, L, maxwidth, cutoff = 10, ...) {

  fn='findPeaks'
  #.msg(sprintf("%s: start.", fn), preLine=FALSE, ...)

  if (cutoff<=0) cutoff=1

  ## if chrnames start with digits (e.g., 1, 2L), then add 'Chr' to prevent errors in definder::regionMatrix
  dc=grep('^[0-9]', names(fullCov))
  changed=FALSE
  if (length(dc)>0) {
    names(fullCov)=paste0('Chr', names(fullCov))
    changed=TRUE
    #.msgin(fn, "Warning: some or all chr names in fullCov start with digits, will add [Chr] to the chr names of the output coverage list!", ...)
  }

  .msgin(fn, sprintf("start regionMatrix (L=%d, cutoff=%d)...", L, cutoff), ...)
  regionMat <- derfinder::regionMatrix(fullCov = fullCov, L = L, cutoff = cutoff)
  ## regionMat: a list, with each chr a Data.Frame(Rle to store coverage at each position)

  .splitPeaks <- function(pos_count, chrom, pos, readlength) {
    chr <- rep(chrom, length(pos_count))
    threshold <- quantile(pos_count, 0.25)
    tempdata <- bumphunter::regionFinder(pos_count, chr = chr, pos = pos, cutoff = threshold, verbose = FALSE)
    tempdata <- subset(tempdata, L > readlength)
    return(tempdata)
  }


  DefPeak <- GenomicRanges::GRanges()
  for (i in names(regionMat)) {
    index <- which(GenomicRanges::width(regionMat[[i]]$regions) > maxwidth)

    if (length(index)>0)
      .msgin(fn, sprintf("chr=%s has %d peaks with width>%dnt, start splitPeaks", i, length(index), maxwidth), ...)

    tempDataFram <- data.frame()
    for (j in index) {
      rawstart <- start(regionMat[[i]]$regions[j]) - 1
      temp <- regionMat[[i]]$bpCoverage[[j]]
      temp$pos <- c(1:dim(temp)[1])
      first.split.peaks <- .splitPeaks(temp$value, i, temp$pos, L)
      pre.split.peaks <- first.split.peaks
      split.peaks <- data.frame()
      iter <- 0
      while (TRUE) {
        if (iter >= 3) {
          break()
        }
        Ind <- which(pre.split.peaks$L > maxwidth)
        cur.split.peaks <- data.frame()
        if (length(Ind) > 0) {
          split.peaks <- rbind(split.peaks, pre.split.peaks[-Ind, ])
          for (k in Ind) {
            pos.temp <- temp[pre.split.peaks$start[k]:pre.split.peaks$end[k], ]
            cur.split.peaks <- .splitPeaks(pos.temp$value, i, pos.temp$pos, L)
          }
        } else {
          split.peaks <- rbind(split.peaks, pre.split.peaks)
          break()
        }
        pre.split.peaks <- cur.split.peaks
        iter <- iter + 1

      }
      split.peaks$start <- split.peaks$start + rawstart
      split.peaks$end <- split.peaks$end + rawstart
      tempDataFram <- rbind(tempDataFram, split.peaks)
    }

    if (length(index)>0) {
      suppressWarnings( DefPeak <- c(DefPeak, regionMat[[i]]$regions[-index]) ) #all+splitedpeaks
    }  else {
      suppressWarnings( DefPeak <- c(DefPeak, regionMat[[i]]$regions) ) #all peaks
    }

    if (nrow(tempDataFram)>0) {
      suppressWarnings(DefPeak <- c(DefPeak, regioneR::toGRanges(tempDataFram)))
      .msgin(fn, sprintf("chr=%s split %d wide peaks to %d small peaks", i, length(index), nrow(tempDataFram)), ...)
    }
    .msgin(fn, sprintf("chr=%s find %d peaks", i, length(DefPeak)), ...)
  }

  #DefPeak <- DefPeak[-which(DefPeak$area == GenomicRanges::width(DefPeak))]
  DefPeak <- DefPeak[GenomicRanges::width(DefPeak) > L]

  if (length(DefPeak)==0) {
    .msgin(fn, sprintf("Warning: 0 peak with width>%dnt (Read Length), please check the readlength (L) parameter!", i, L), ...)
  } else {
    .msgin(fn, sprintf("%d peaks with width>%dnt (Read Length) found on all chrs", length(DefPeak), L), ...)
  }

  strand(DefPeak) <- strand
  if (length(DefPeak)==0) {
    DefPeak=as.data.frame(DefPeak)
  } else {
    DefPeak=as.data.frame(DefPeak, row.names = 1:length(DefPeak)) #avoid dup row.names
  }
  DefPeak <- dplyr::mutate(DefPeak, chr = seqnames)
  DefPeak=DefPeak[, c('chr','strand','start','end','value')]

  # strand(DefPeak) <- strand
  # if (strand == "+") {
  #     DefPeak <- as.data.frame(DefPeak, row.names = 1:length(DefPeak))
  #     DefPeak <- dplyr::mutate(DefPeak, chr = seqnames, coord = end)
  # } else {
  #     DefPeak <- as.data.frame(DefPeak, row.names = 1:length(DefPeak))
  #     DefPeak <- dplyr::mutate(DefPeak, chr = seqnames, coord = start)
  # }

  # restore chr names
  if (changed & nrow(DefPeak)>0) {
    #.msgin(fn, "Restore chr names by removing [Chr]", ...)
    DefPeak$chr=gsub('^Chr','',DefPeak$chr)
  }

 # .msg(sprintf("%s: finish.", fn), ...)
  return(DefPeak)
}

write.bed <- function(.x, f) {
    write.table(x = .x, file = f, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}


#' Perform peak calling on a BAM file with single strand
#'
#' findPeaksByStrand calls peaks on a BAM file with single strand. This function combines two functions in previous scAPAtrap, loadBpCoverages and findPeaks.
#'
#' @inheritParams findPeaks
#' @param bamFile a BAM file with index (.bai)
#' @param chrs Chromosome information. If not provided, then will try get chrs from the bamFile.
#' @param ofile output file with each line is one peak
#' @return If ofile=NULL, return a data frame with many columns but four fixed columns chr/strand/start/end/value, with each row being a peak.
#' Otherwise, output to file and return ofile name.
#' @examples
#' \dontrun{
#' forwardPeaks <-findPeaksByStrand(bamFile, chrs=NULL, strand='+', L=98, maxwidth=1000, cutoff=10)
#' reversePeaks <-findPeaksByStrand(bamFile, chrs=NULL, strand='-', L=98, maxwidth=1000, cutoff=10)
#' }
#' @export
findPeaksByStrand <- function(bamFile, chrs=NULL, strand, L, maxwidth, cutoff = 10, ofile=NULL, ...) {

  fn='findPeaksByStrand'
  .msg(sprintf("%s (%s): start.", fn, strand), preLine=FALSE, ...)

  if (cutoff<=0) cutoff=1

  if (length(chrs)==0) {
    chrs=.getBAMchrs(bamFile)
    .msgin(fn, string=sprintf("Total %d chrs (e.g., %s)", length(chrs), .shortStr(chrs)), ...)
  }

  if (length(chrs)==0) {
    .msgin(fn, string="Get chr error in bamFile!", ...)
    stop("findPeaksByStrand error, please check logf to see Logs if set logf!")
  }

  .msgin(fn, sprintf("loadBpCoverages on %d chrs: %s", length(chrs), .shortStr(chrs)), ...)
  fullCov <- derfinder::fullCoverage(files = bamFile, chrs = chrs, ...)
  ## fullCov: a list, with each chr a Data.Frame(Rle to store coverage at each position)

  DefPeak=findPeaks(fullCov=fullCov, strand=strand, L=L, maxwidth=maxwidth, cutoff = cutoff, ...)

  if (!is.null(ofile)) {
    write.table(DefPeak, file=ofile, col.names=TRUE, row.names = FALSE, sep="\t", quote=FALSE)
    .addTrapFiles(ofile, paste0('*findPeaksByStrand: peak (', strand, ')'))
    .msg(sprintf("%s: finish.", fn), ...)
    return(ofile)
  } else {
    .msg(sprintf("%s (%s): finish.", fn, strand), ...)
    return(DefPeak)
  }
}

#----------------------------------- countPeaks ----------------------------------------------
#' Generate a peak annotation file
#'
#' @description Generate a peak annotation file (PeakID/chr/start/end/Strand) by combining forward peaks and reverse peaks generated by findPeaks
#'
#' @param forwardPeaks A peak filename or a peak data.frame of forward strand, which is the output of \code{\link{findPeaksByStrand}} or \code{\link{findPeaks}}, where strand = +.
#' @param reversePeaks Save as forwardPeaks but on the reverse strand. Can provide forwardPeaks or reversePeaks, or both.
#' @param outputdir Output file directory.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{If `TRUE` basic status updates will be printed along the way.}
#' \item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
#' }
#' @examples
#' \dontrun{
#' forwardPeaks <-findPeaks(fullcovF,'+', 98, 1000, cutoff = 1)
#' reversePeaks <-findPeaks(fullcovR,'-',98, 1000, cutoff = 1)
#' generateSAF(forwardPeaks, reversePeaks,'./data')
#' }
#' @return Full path of the peak list (<outputdir>/peaks.saf). This file contains five columns without header: peakID/chr/start/end/strand.
#' @export
generateSAF <- function(forwardPeaks, reversePeaks, outputdir, ...) {

  fn='generateSAF'
  .msg(sprintf("%s: start.", fn), preLine=FALSE, ...)

  if (missing(x = forwardPeaks) && missing(x = reversePeaks)) {
    stop("Must provide at least forwardPeaks or reversePeaks (data.frame or filename)!")
  }

  both=FALSE
  if (!missing(x = forwardPeaks) && !missing(x = reversePeaks)) both=TRUE

  .loadFRpeaks<-function(frfile) {
    if (is.character(frfile)) {
      if (!file.exists(frfile))
        stop(frfile, 'not exists!')
      frfile <- utils::read.delim(frfile, header = TRUE) # chr/strand/coord/start/end... not known
    }
    return(frfile)
  }

  if (!is.null(forwardPeaks))
    forwardPeaks=.loadFRpeaks(forwardPeaks)
  if (!is.null(reversePeaks))
    reversePeaks=.loadFRpeaks(reversePeaks)
  if (both) peaks <- rbind(.loadFRpeaks(forwardPeaks), .loadFRpeaks(reversePeaks))

  peaks$PeakID <- paste0("peak", "_", 1:length(peaks$chr))
  peaks.saf <- peaks[, c("PeakID", "chr", "start", "end", "strand")]

  if (!file.exists(outputdir)) {
      dir.create(outputdir)
  }

  output <- paste0(outputdir, "/", "peaks.saf")
  write.bed(peaks.saf, output)
  .addTrapFiles(output, "generateSAF peaks.saf")

  .msg(sprintf("%s: finish.", fn), ...)
  return(output)
}

.mergePA <- function(bam, adjBp = 10) {
    peak <- GenomicRanges::GRanges(seq = S4Vectors::Rle(bam$chr), ranges = IRanges::IRanges(bam$coord, bam$coord + adjBp), strand = S4Vectors::Rle(bam$strand))

    peak <- GenomicRanges::reduce(peak)
    peak <- as.data.frame(peak)
    peak$end <- peak$end - adjBp
    peak$width <- peak$end - peak$start + 1

    peak$coord <- 0
    peak$coord[peak$strand == "+"] <- peak$end[peak$strand == "+"]
    peak$coord[peak$strand == "-"] <- peak$start[peak$strand == "-"]
    names(peak) <- c("chr", names(peak)[-1])
    return(peak)
}


.computePAcoord <- function(peak) {
    peak$coord <- 0
    peak$coord[peak$strand == "+"] <- peak$end[peak$strand == "+"]
    peak$coord[peak$strand == "-"] <- peak$start[peak$strand == "-"]

    return(peak)
}


#' Generate a BAM file for called peaks with featureCounts and samtools
#'
#' @description Generate a BAM file for called peaks with featureCounts and samtools, which can be used by umi_tools to quantify peaks.
#' @param featureCounts.path The path of the featureCounts tool.
#' @param samtools.path The path of the samtools.
#' @param input BAM file.
#' @param peakfile The path of the peaks.saf file generated by \code{\link{generateSAF}}.
#' @param thread Number of CPU threads, default is 12.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{If `TRUE` basic status updates will be printed along the way.}
#' \item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
#' }
#' @examples
#' \dontrun{
#' input <- './data/demo.bam'
#' peakfile <- './data/peaks.saf'
#' generateFinalBam(featureCounts.path, samtools.path, input, peakfile, 24)
#' }
#' @return Generate final.bam file (final.bam) in the folder of the <input> dir.
#' @export
generateFinalBam <- function(featureCounts.path, samtools.path, input, peakfile, thread=12, ...) {

  fn='generateFinalBam'
  .msg(sprintf("%s: start.", fn), preLine=FALSE, ...)

  output.dir <- dirname(input)
  log <- paste0(output.dir, "/peak_assigned")
  command <- paste0(featureCounts.path, " -a ", peakfile, " -F SAF -t exon -s 1 -M --largestOverlap -o ", log, " -R BAM -T ",
      thread, " ", input)
  ofile=paste0(input, ".featureCounts.bam")
  .runCommand(fn, command, ofile=ofile, ocheck=TRUE, ...)

  .addTrapFiles(log, '*generateFinalBam: featureCounts log')
  .addTrapFiles(ofile, '*generateFinalBam: featureCounts BAM')

  input <- paste0(input, ".featureCounts.bam")
  output <- paste0(dirname(input), "/", "final.bam")
  command <- paste0(samtools.path, " sort -@ ", thread, " ", input, " -o ", output, " && ", samtools.path, " index -@ ",
      thread, " ", output)

  .runCommand(fn, command, ofile=output, ocheck=TRUE, ...)
  .addTrapFiles(output, '*generateFinalBam: samtools sort')

  .msg(sprintf("%s: finish.", fn), ...)

  return(output)
}

#' Calculate the expression level of each peak in each cell with umi_tools
#'
#' @param umitools.path The path of the umi_tools.
#' @param input The final.bam file is generated by \code{\link{generateFinalBam}}.
#' @param outputdir Output file directory.
#' @param TenX Logical value, TRUE for 10X data or BAM resulted from STARsolo.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{If `TRUE` basic status updates will be printed along the way.}
#' \item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
#' \item{notRun }{Default is FALSE. If `TRUE`, the Shell commands inside this function will not run but output the commad line.}
#' }
#' @examples
#' umitools.path <- '/home/aa/miniconda2/envs/umi_tools/bin/umi_tools'
#' input <- './data/final.bam'
#' outputdir <- './data'
#' \dontrun{
#' countPeaks(umitools.path, input, outputdir, TenX=TRUE)
#' }
#' countPeaks(umitools.path, input, outputdir, TenX=TRUE, notRun=TRUE)
#' @return Full path of the <counts.tsv.gz> file in the outputdir directory
#' @export
countPeaks <- function(umitools.path, input, outputdir, TenX = TRUE, ...) {

  fn='countPeaks'
  .msg(sprintf("%s: start.", fn), preLine=FALSE, ...)

  output <- paste0(outputdir, "/counts.tsv.gz")
  if (TenX) {
    command <- paste0(umitools.path, " count --per-gene --gene-tag=XT --assigned-status-tag=XS --method=unique --extract-umi-method=tag --umi-tag=UB --cell-tag=CB --per-cell -I ",
        input, " -S ", output)
  } else {
    command <- paste0(umitools.path, " count --per-gene --gene-tag=XT --assigned-status-tag=XS --method=unique --per-cell -I ",
        input, " -S ", output)
  }
  .runCommand(fn, command, ofile=output, ocheck=TRUE, ...)
  .addTrapFiles(output, 'countPeaks: umitools counts')
  .msg(sprintf("%s: finish.", fn), ...)
  return(output)
}

.findIntersectBypeak <- function(qrypeak, sbjpeak, d = 100) {
  gr1 <- with(qrypeak, GRanges(seqnames = chr, ranges = IRanges::IRanges(start = coord, end = coord), strand = strand))
  gr2 <- with(sbjpeak, GRanges(seqnames = chr, ranges = IRanges::IRanges(start = coord, end = coord), strand = strand))
  ov = GenomicRanges::findOverlaps(gr1, gr2, maxgap = d-1, minoverlap = 0L, type = c("any"), select = "all", ignore.strand = FALSE)
  return(qrypeak[unique(ov@from), ])
}

#----------------------------------- output scAPAtrapData ----------------------------------------------
# validate.arg(NULL, c('a','b','c', 'no')) --> a
# validate.arg('x', c('a','b','c', 'no')) --> error
# validate.arg(NULL, c('a','b','c', 'no'), null2default=FALSE)  --> error
# validate.arg(NULL, c('a','b','c', 'no', NULL), null2default=FALSE)  --> error
# validate.arg('a', c('ca','ab','c', 'no')) --> ab
# validate.arg('A', c('ca','ab','c', 'no')) --> ab
validate.arg <- function(arg, choices, several.ok = FALSE, null2default = TRUE, lc = TRUE) {
  if (is.null(arg)) {
    if (null2default) arg=choices[1]
  }
  if (lc) arg=tolower(arg)
  arg = match.arg(arg, choices, several.ok)
  return(arg)
}

## load a count table or read countsfile, and check
## return 3-column count table (gene, cell, count)
.loadCounts<-function(countsfile) {
  if (is.character(countsfile)) {
    if (!file.exists(countsfile))
      stop(countsfile, 'not exists!')
    countsfile <- utils::read.delim(countsfile, header = TRUE) # gene (peakID), cell, count
  }

  if (ncol(countsfile)!=3) stop("countsfile should be of 3 columns: gene/cell/count")
  if (!all(c('gene','cell','count') %in% colnames(countsfile))) stop("countsfile should be of 3 columns: gene/cell/count")
  return(countsfile)
}

.loadPeaks<-function(peaksfile) {
  if (is.character(peaksfile)) {
    if (!file.exists(peaksfile))
      stop(peaksfile, 'not exists!')
    peaksfile <- utils::read.delim(peaksfile, header = FALSE)
    names(peaksfile) <- c("peakID", "chr", "start", "end", "strand")
    rownames(peaksfile)=peaksfile$peakID
  }

  if (ncol(peaksfile)!=5) stop("peaks should be of 5 columns: peakID/chr/start/end/strand")
  if (!all(c('peakID','chr','start','end','strand') %in% colnames(peaksfile))) stop("peaks should be of 5 columns: peakID/chr/start/end/strand")
  rownames(peaksfile)=peaksfile$peakID
  return(peaksfile)
}


.loadTails<-function(tailsfile) {
  if (is.character(tailsfile)) {
    if (!file.exists(tailsfile))
      stop(tailsfile, 'not exists!')
    tailsfile <- utils::read.delim(tailsfile, header = TRUE) # chr/strand/coord/count
  }

  if (ncol(tailsfile)!=4) stop("tailsfile should be of 4 columns: chr/strand/coord/count")
  if (!all(c('chr','strand','coord','count') %in% colnames(tailsfile))) stop("tailsfile should be of 4 columns: chr/strand/coord/count")
  return(tailsfile)
}


#' Reduce peaks number in a countsfile and/or peaksfile
#'
#' @description Reduce peaks in a countsfile or counts table and also remove same peaks in peaksfile (if provided), by min.cells/max.cells, and min.counts/max.counts.
#' This function is useful to retrieve highly expressed, lowly expressed peaks or moderately expressed peaks.
#'
#' @param countsfile The decompressed file path of counts.tsv.gz generated by \code{\link{countPeaks}}, or the count table with three columns.
#' @param peaksfile peaksfile or peak table with five columns. If not NULL, then filter peaksfile after filtering countsfile.
#' @param min.cells retain peaks expressed in >= min.cells, the default value is 10.
#' @param min.count retain peaks with read count >= min.count, the default value is 10.
#' @param max.cells retain peaks expressed in < max.cells, the default value is NULL (unlimited). This is used to filter peaks with less expression.
#' @param max.count retain peaks with read count < max.count, the default value is NULL (unlimited). This is used to filter peaks with less expression.
#' @param suffix applicable when countsfile and peaksfile are both provided. Then counts and peaks will be output to <countsfile>.reduced; <peaksfile>.reduced.
#' @param toSparse to output a sparseMatrix (gene-cell) or keep the triplet table as input.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{If `TRUE` basic status updates will be printed along the way.}
#' \item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
#' }
#' @return A data.frame (toSparse=FALSE), or a sparse Matrix (toSparse=TRUE) of counts,
#' or a filename list with (countsfile, peaksfile) (if peaksfile is not NULL).
#' @examples
#' \dontrun{
#' countsfile='../dataFly/APA.tails.no/counts.tsv.gz'
#' peaksfile='../dataFly/APA.tails.no/peaks-notails.saf'
#' ## only filter countsfile or counts-table, return a df
#' reducePeaks(countsfile, min.cells = 10, min.count = 10, toSparse=FALSE)
#'
#' ## retain large peaks, and output both counts and peaks, save to .reduced file (>=10 & >=50)
#' reducePeaks(countsfile=countsfile, peaksfile=peaksfile, min.cells = 10, min.count = 50)
#'
#' ## retain low-expressed peaks, and output both counts and peaks, save to .reduced file (<=9 & <=49)
#' reducePeaks(countsfile=countsfile, peaksfile=peaksfile, max.cells = 9, max.count = 49,
#'             min.cells=NULL, min.count=NULL, suffix='.small')
#'
#' <=9
#' reducePeaks(countsfile=countsfile, peaksfile=peaksfile, max.cells = 9,  max.count=NULL,
#'            min.cells=NULL, min.count=NULL, suffix='.smallcells')
#'
#' <=49
#' reducePeaks(countsfile=countsfile, peaksfile=peaksfile, max.cells = NULL,  max.count=49,
#'             min.cells=NULL, min.count=NULL, suffix='.smallcounts')
#'
#' smallpeaks=.loadPeaks('../dataFly/APA.tails.no/peaks-notails.saf.small')
#' largepeaks=.loadPeaks('../dataFly/APA.tails.no/peaks-notails.saf.reduced')
#' smallpeaks1=.loadPeaks('../dataFly/APA.tails.no/peaks-notails.saf.smallcells')
#' smallpeaks2=.loadPeaks('../dataFly/APA.tails.no/peaks-notails.saf.smallcounts')
#' fullpeaks=.loadPeaks(peaksfile)
#' nrow(smallpeaks); nrow(largepeaks); nrow(fullpeaks)
#' smallset2=unique(rbind(smallpeaks1, smallpeaks2))
#' nrow(smallset2) + nrow(largepeaks);  nrow(fullpeaks)
#' ## should be the same, but if not the same, may be some peakIDs in fullpeaks are not in the counts table
#'
#' }
#' @export
reducePeaks <- function(countsfile, peaksfile=NULL,
                        min.cells = 10, min.count = 10,
                        max.cells = NULL, max.count = NULL,
                        suffix='.reduced', toSparse=FALSE, ...) {

  if (!is.null(peaksfile)) {
    if (!is.character(countsfile) & !is.character(peaksfile)) stop("reducePeaks: countsfile and peaksfile should be filename when both are provided\n")
    if (!file.exists(countsfile) | !file.exists(peaksfile)) stop("reducePeaks: countsfile or peaksfile not exist!")
    if (suffix=='') stop("reducePeaks: countsfile and peaksfile both are provided, suffix (e.g., 'reduced' to the filename) should be provided\n")
    toSparse=FALSE
  }

  fn='reducePeaks'
  counts=.loadCounts(countsfile)
  N1=nrow(counts)

  counts = transform(counts,
                     gene = factor(gene),
                     cell = factor(cell))
  scounts = Matrix::sparseMatrix(as.integer(counts$gene), as.integer(counts$cell), x = counts$count)
  colnames(scounts) = levels(counts$cell)
  rownames(scounts) = levels(counts$gene)
  invisible(gc())

  N=nrow(scounts)

  nc=Matrix::rowSums(scounts) #counts
  nv=Matrix::rowSums(scounts>0) #cells

  nids=list(min.count=NA, min.cells=NA, max.count=NA, max.cells=NA)
  if (!is.null(min.count)) nids[[1]]=which(nc>=min.count)
  if (!is.null(min.cells)) nids[[2]]=which(nv>=min.cells)
  if (!is.null(max.count)) nids[[3]]=which(nc<=max.count)
  if (!is.null(max.cells)) nids[[4]]=which(nv<=max.cells)

  for (i in 1:length(nids)) {
    if (!is.na(nids[[i]][1])) {
      .msgin(fn, sprintf("%d peaks meet %s (%d).", length(nids[[i]]), names(nids)[i], get(names(nids)[i])), ...)
    }
  }

  ## get intersection of all conditions
  nids[is.na(nids)]=NULL
  allids=nids[[1]]
  if (length(nids)>=2) {
    for (i in 2:length(nids)) {
      allids=intersect(allids, nids[[i]])
    }
  }

  scounts=scounts[allids, ]
  .msgin(fn, sprintf("%d peaks (total %d) meet combined conditions", length(allids), N), ...)

  if (toSparse) return(scounts)

  peakid=rownames(scounts) ##sparse rowname = peakid

  #convert back to triplet table
  rid=rownames(scounts)
  cid=colnames(scounts)
  scounts <- as.data.frame(summary(scounts))
  colnames(scounts)=c('i','j','count')
  scounts$gene <- rid[scounts$i]
  scounts$cell <-cid[scounts$j]
  scounts=transform(scounts, i=NULL, j=NULL)
  scounts=scounts[, c('gene','cell','count')]

  if (is.null(peaksfile)) return(scounts)

  ## save both peaks and counts
  if (!is.null(peaksfile)) {
    saf=.loadPeaks(peaksfile)
    ipeak=base::intersect(saf$peakID, peakid) # common peaks
    if (length(ipeak)!=length(peakid)) {
      .msgin(fn, sprintf("Warning: peakID in counts-table are not all in peaksfile, will use common peaks (%d) to filter both counts (%d) and peaks (%d)!",
                         length(ipeak), length(peakid), nrow(saf)), ...)
    }
    if (!(all(ipeak %in% saf$peakID))) {
      .msgin(fn, sprintf("Warning: peakID in peaksfile are not all in counts-table, will use common peaks (%d) to filter both counts (%d) and peaks (%d)!",
                         length(ipeak), length(peakid), nrow(saf)), ...)
    }
    saf=saf[ipeak, ]
    scounts=scounts[scounts$gene %in% ipeak, , drop=FALSE]
    write.table(saf, file=paste0(peaksfile, suffix), col.names = FALSE, row.names = FALSE, sep="\t", quote=FALSE)
    write.table(scounts, file=paste0(countsfile, suffix), col.names = TRUE, row.names = FALSE, sep="\t", quote=FALSE)

    .addTrapFiles(paste0(peaksfile, suffix), 'reducePeaks: reduced peaksfile')
    .addTrapFiles(paste0(countsfile, suffix), 'reducePeaks: reduced countsfile')

    .msgin(fn, sprintf("save reduced peaksfile to <%s> and countsfile to <%s>",
                       basename(paste0(peaksfile, suffix)), basename(paste0(countsfile, suffix)) ), ...)
    return(list(countsfile=paste0(countsfile, suffix), peaksfile=paste0(peaksfile, suffix)))
  }
}

#' Generate single cell expression matrix
#'
#' @description Given peaks' count and annotation file, this function can filter peaks with low expression levels and generate different types of objects,
#'              including data.frame (df), list with data.frame peaks.meta and Sparse Matrix peaks.count (list), SeuratObject (Seurat),  SingleCellExperiment (sce), movAPA's PACdataset (PAC).
#'
#' @param countsfile The decompressed file path of counts.tsv.gz generated by \code{\link{countPeaks}} or a count table (gene/cell/count).
#' @param peaksfile The path of the peaks.saf file generated by \code{\link{generateSAF}}.
#' @param barcode Vector of the barcode list for filtering cells.
#' @param tails The result produced by \code{\link{findTails}} or \code{\link{findTailsByPeaks}}.
#' If provided, then filter tailed-peaks within d nt.
#' @param d Max distance from the tail and peak to determine whether a peak is real (if it is within d-nt of a tail).
#' @param min.cells peak is expressed in at least a few cells,the default value is 10.
#' @param min.count peak minimum expression, the default value is 10.
#' @param ofile if not NULL then output the an variable called scAPAtrapData to ofile, and return ofile.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{If `TRUE` basic status updates will be printed along the way.}
#' \item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
#' }
#' @return A list containing peaks.meta and peaks.count (ofile=NULL) or the ofile name (ofile!=NULL).
#' The ofile is an RDA file which contains a variable called scAPAtrapData after load(ofile). The scAPAtrapData is a list.
#' @examples
#' \dontrun{
#' countsfile <- countPeaks(umitools.path,'./data/final.bam','./data',TenX=TRUE)
#' peakfile <- generateSAF(forwardPeaks,reversePeaks,'./data')
#' tails <- findTails(bamfile = './data/demo.bam')
#' barcode <- utils::read.delim2('./data/barcodes.tsv',header = FALSE)
#' barcode <- gsub('-[0-9]','',barcode$V1)
#' scExpMa <- generatescExpMa(countsfile, peaksfile, barcode, tails,
#'                           d=100, min.cells = 10, min.count = 10)
#' scExpMa2 <- eneratescExpMa(countsfile, peaksfile, barcode,
#'                min.cells = 10, min.count = 100,  verbose=TRUE, logf='xx.log')
#' }
#' @section Warning: The chromosome name in tails and peaksfile(peaks.saf) must be consistent (e.g., both be 1 or be Chr1).
#' @export
generatescExpMa <- function(countsfile, peaksfile, barcode=NULL,
                            tails=NULL, d=100,
                            min.cells = 10, min.count = 10,
                            ofile=NULL, ...) {

  fn='generatescExpMa'

  .msg(sprintf("%s: start.", fn), preLine=FALSE, ...)

  counts=.loadCounts(countsfile)
  saf=.loadPeaks(peaksfile)
  invisible(gc())

  # filter non-barcode cells in counts
  if (!is.null(barcode)) {
    cs1=length(unique(counts$cell))
    counts <- subset(counts, cell %in% barcode)
    cs2=length(unique(counts$cell))
    if (cs1!=cs2) {
      .msgin(fn, sprintf("%d cells from %d total cells are in barcode (#%d).", cs2, cs1, length(barcode)), ...)
    }
  }

  # add coord
  saf <- .computePAcoord(saf)

  # use tail-positions to filter peaks within d-nt
  if (!is.null(tails)) {
    tails=.loadTails(tails)
    .msgin(fn, sprintf("Using %d tails to filter tailed-peaks within %d nt.", nrow(tails), d), ...)
    n1=nrow(saf)
    saf <- .findIntersectBypeak(saf, tails, d=d)
    n2=nrow(saf)
    .msgin(fn, sprintf("%d non-tailed peaks are removed.",  n1-n2), ...)
    invisible(gc())

    ## filter counts
    n1=nrow(counts)
    counts=counts[counts$gene %in% saf$peakID, ]
    n2=nrow(counts)
    .msgin(fn, sprintf("%d non-tailed peak-count rows are removed", n1-n2), ...)
  }

  ## slow BIG MEM!!
  ##scExpMa <- reshape2::dcast(counts, gene ~ cell, fun.aggregate=sum, value.var = "count")

  # get gene-cell matrix very fast with Sparse Matrix
  #https://www.r-bloggers.com/2016/01/casting-a-wide-and-sparse-matrix-in-r/
  .msgin(fn, 'get peak-cell count matrix.', ...)

  scExpMa=reducePeaks(counts, min.cells = min.cells, min.count = min.count, toSparse=TRUE, ...)

  #make counts and meta data in order, but in different matrices
  .msgin(fn, "link peaks.meta and peaks.count.", ...)
  n1=nrow(saf)
  n2=nrow(scExpMa)
  cmpeaks=base::intersect(saf$peakID, rownames(scExpMa))
  n3=length(cmpeaks)
  .msgin(fn, sprintf("peaks.meta=%d peaks; peaks.counts=%d peaks; common=%d peaks.", n1, n2, n3), ...)

  rownames(saf)=saf$peakID
  saf=saf[cmpeaks, , drop=FALSE]
  scExpMa=scExpMa[cmpeaks, ,drop=FALSE]
  if (!identical(saf$peakID, rownames(scExpMa))) stop('Peaks in peaks.meta are not the same as the peak-cell matrix!\n')

  .msgin(fn, sprintf('output a list with %d peaks.meta (data.frame) and peaks.count (Sparse Matrix).', nrow(saf)), ...)
  if (!is.null(ofile)) {
    scAPAtrapData=list(peaks.meta=saf, peaks.count=scExpMa)
    save(scAPAtrapData, file=ofile)
    .addTrapFiles(ofile, 'generatescExpMa: scAPAtrapData')
    .msg(sprintf("%s: finish.", fn), ...)
    return(ofile)
  } else {
    .msg(sprintf("%s: finish.", fn), ...)
    return(list(peaks.meta=saf, peaks.count=scExpMa))
  }
}

#' Convert the scExpMa list generated by scAPAtrap to other format
#'
#' @description Given peaks' count and annotation file, this function can filter peaks with low expression levels and generate different types of objects,
#'              including data.frame (df), list with data.frame peaks.meta and Sparse Matrix peaks.count (list), SeuratObject (Seurat),  SingleCellExperiment (sce), movAPA's PACdataset (PAC).
#'
#' @param traplist A list with data.frame peaks.meta and Sparse Matrix peaks.count generated by generatescExpMa.
#' @param oformat A character string of df/seurat/sce/pac to generate different types of objects:
#'                 data.frame (df),  SeuratObject (Seurat),  SingleCellExperiment (sce), movAPA's PACdataset (PAC).
#'                 Generally, the output includes peaks.count (SparseMatrix) and peaks.meta (peakID/chr/start/end/strand/coord).
#'                 oformat=df may be problematic when the peak-cell matrix is very large. If oformat=df failed, please set oformat to other strings to represent the peak-cell counts matrix as Sparse Matrix.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{If `TRUE` basic status updates will be printed along the way.}
#' \item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
#' }
#' @examples
#' \dontrun{
#' convertAPAtrapData(from=peakslist, to='pac')
#' }
#' @export
convertAPAtrapData<-function(traplist, oformat, ...) {

  fn='convertAPAtrapRes'

  oformat=validate.arg(oformat, c('df', 'seurat', 'pac', 'sce'))

  if (is.character(traplist)) {
    if (!file.exists(traplist)) stop('convertAPAtrapData: traplist is a filename but not exists!')
    e=new.env()
    v=load(traplist, envir = e)
    traplist=get(v, envir=e)
    rm(list=v, envir=e)
    .msgin(fn, paste0('traplist is an .rda file, loaded scAPAtrapData.'), ...)
  }

  if (!is.list(traplist)) stop("convertAPAtrapRes: traplist is not a list!")
  if (!(all(c('peaks.meta','peaks.count') %in% names(traplist)))) stop("convertAPAtrapRes: peaks.meta or peaks.count not in traplist!")

  scExpMa=traplist$peaks.count
  saf=traplist$peaks.meta

  if (oformat=='df') {
    msg="Warning: oformat=df may be problematic when the peak-cell matrix is very large.\nIf failed, please set pformat=list/seurat to use sparse matrix for storing the peak-count table!"
    .msgin(fn, paste0('oformat=df, output a data frame with peaks.meta and peaks.count together.', msg), ...)
    scExpMa=as.data.frame(as.matrix(scExpMa))
    invisible(gc())
    saf=cbind(saf, scExpMa)
    return(saf)
  }

  if (oformat=='seurat') {
    .msgin(fn, 'oformat=seurat, output a SeuratObject (assay=counts).', ...)
    so<-Seurat::CreateSeuratObject(scExpMa, assay = "counts")

    # make feature names like "peak-XXXX" , same as the assay in Seurat obj
    rownames(saf)=gsub('_', '-', rownames(saf))
    so[['counts']]=Seurat::AddMetaData(so[['counts']], saf)

    # Access cell-level meta-data
    #head(so[[]])
    # Access feature-level meta-data
    #head(so[["counts"]][[]])
    # Access counts table
    #so[["counts"]][1:10, 1:5]
    return(so)
  }

  if (oformat=='sce') {
    .msgin(fn, 'oformat=sce, output SingleCellExperiment (assay=counts).', ...)
    sce=SingleCellExperiment::SingleCellExperiment(assays=list(counts=scExpMa))
    rowData(sce)=saf
    return(sce)
    #head(rowData(sce))
    #head(counts(sce)[1:10, 1:10])
    #head(assay(sce, "counts")[1:10, 1:10])
  }

  if (oformat=='pac') {
    .msgin(fn, 'oformat=pac, output movAPA::PACdataset.', ...)
    pacds=movAPA::createPACdataset(counts=scExpMa, anno=saf)
    return(pacds)
  }

}


############# scAPAtrap wrapper #############
## global variable TRAP.PARAMS to store scAPAtrap parameters.
## default values
myenv$TRAP.PARAMS=list(
  TenX=TRUE, #whether 10X or STARsolo BAM
  barcode=NULL, #barcode vectors
  chrs=NULL, #chromosome names, NULL to retrieve from BAM
  maxwidth=1000,  #max pA peak width
  readlength=49,  #R2 read length
  cov.cutoff=10, #min coverage for peak calling
  min.cells=10, #filter peaks present in >=10 cells
  min.count=10, #filter peaks with >=10 reads
  tails.search='no', #search tails
  d=100, #min distance from tail to peak
  findUniqueMap=TRUE, #whether to run findUniqueMap
  findUniqueMap.sort=TRUE, #whether to sort BAM in findUniqueMap
  findUniqueMap.index=TRUE, #whether to index BAM in findUniqueMap
  thread=24 #thread CPU
)

#' Get parameters to run scAPAtrap
#'
#' @details
#' TRAP.PARAMS returns the global variable of list type to store parameters for running the wrapper function \code{\link{scAPAtrap}}.
#' \itemize{
#' \item{TenX}: Whether the inputBam is 10X-format (e.g., BAM file from CellRanger or STARsolo).
#' For non-10X protocols like CEL-seq, when STARsolo is used, the BAM file will also contain cell barcode tag (CB) and UMI barcode tag (UB), then can set TenX=TRUE.
#' \item{barcode}: string vector to define barcodes, default is NULL (not provided).
#' \item{chrs}: a string vector of chromosome names (e.g., chr1, chr2, chr3...), default is NULL.
#' This parameter is used in \code{\link{findPeaksByStrand}}.
#' If chrs=NULL and run one-step \code{\link{scAPAtrap}}, then will try to get all chrs from inputBAM.
#' However, when running \code{\link{scAPAtrap}}, if chrs=NULL and there are more than 100 chrs from inputBAM,
#' users should provide trap.params$chrs explicitly to run scAPAtrap just in case there are too many unuseful chrs to run.
#' \item{maxwidth}: max width (nt) to define a valid peak, default is 1000.
#' Peaks wider than `maxwidth` will be automatelly split until meet requirement.
#' This parameter is used in \code{\link{findPeaksByStrand}}.
#' \item{readlength}: R2 read length, default is 49. `readlength` determines the min peak width;
#' peaks wider than `readlength` will be retained.
#' This parameter is used in \code{\link{findPeaksByStrand}}.
#' \item{cov.cutoff}: min read coverage for peak calling, default is 10.
#' Peak regions with coverage greater than `cov.cutoff`will be considered as candidate peaks.
#' This parameter is used in \code{\link{findPeaksByStrand}}.
#' \item{min.cells}: number of cells that a peak is expressed (i.e., with read count over `min.count`), default is 10.
#' `min.cells` is used to filter peaks present in at least `min.cells` cells.
#' This parameter is used in \code{\link{reducePeaks}} and \code{\link{generatescExpMa}}.
#' \item{min.count}: number of read counts a peak has, default is 10.
#' `min.count` is used to filter peaks with at least `min.count' reads.
#' This parameter is used in \code{\link{reducePeaks}} and \code{\link{generatescExpMa}}.
#' \item{tails.search}: the strategy to search A tails, can be set `peaks`, `genome`, `no` (default). This parameter is used in the wrapper function \code{\link{scAPAtrap}}.
#' \describe{
#' \item{no }{Not search tails, which means the final peaks are not linked with any A-tails.}
#' \item{genome }{Search tails genome-wide, which may be suitable for small genome. Otherwise it will take very long time or MEM.}
#' \item{peaks }{Search tails within the peak regions by `d` nt. This way may be suitable for big genome but with relatively small number of peaks.}
#' }
#' \item{d}: min distance from a peak to any tail to call a peak as real PA, default is 100.
#' Only applicable when tails.search is peaks or genome. This parameter is used in \code{\link{findTailsByPeaks}} and \code{\link{generatescExpMa}}.
#' If provided, then the final scAPAtrapData (the PA list) will remove cells not in `barcode.`
#' \item{findUniqueMap}: default is TRUE, means to run \code{\link{findUniqueMap}} to get unique mappings, and/or sort/index the inputBam.
#' \item{findUniqueMap.sort}: default is TRUE, to sort inputBam. Only applicable when findUniqueMap=TRUE.
#' \item{findUniqueMap.index}: default is TRUE, to index inputBam. Only applicable when findUniqueMap=TRUE.
#' \item{thread}: Number of CPU threads, default is 24.
#' }
#' @return A list of parameters with default values for running \code{\link{scAPAtrap}}.
#' @examples
#' ## see TRAP.PARAMS
#' TRAP.PARAMS()
#' ## set up new trap.params
#' trap.params=TRAP.PARAMS()
#' trap.params$TenX=FALSE
#' ## set up from scratch
#' trap.params=list(thread=24,
#'                  TenX=TRUE,
#'                  chrs=c('2L','2R','3L','3R','4','X','Y'),
#'                  maxwidth=1000,
#'                  readlength=49,
#'                  cov.cutoff=10,
#'                  min.cells=10,
#'                  min.count=10,
#'                  tails.search='peaks',
#'                  d=100,
#'                  barcode=NULL,
#'                  findUniqueMap=TRUE,
#'                  findUniqueMap.sort=TRUE,
#'                  findUniqueMap.index=TRUE )
#' @export
TRAP.PARAMS<-function() {
  return(myenv$TRAP.PARAMS)
}

#' Set the parameters for running scAPAtrap
#'
#' setTrapParams sets parameters for running scAPAtrap, by setting default values using \code{\link{TRAP.PARAMS}} and checking the validity.
#'
#' @param trap.params A list similar to \code{\link{TRAP.PARAMS}}.
#' @param print Whether to print the list (default is TRUE).
#' @param check Whether to check the validity of values (default is TRUE).
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{If `TRUE` basic status updates will be printed along the way.}
#' \item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
#' }
#' @return A list like code{\link{TRAP.PARAMS}}.
#' @examples
#' ## get default
#' trap.params=setTrapParams()
#' ## change some parameters
#' trap.params=setTrapParams(list(chrs=c('2L','2R','3L','3R','4','X','Y'), TenX=FALSE))
#' ## will ignore non-valid parameter names
#' trap.params=setTrapParams(list(xx=FALSE))
#' ## add barcodes
#' trap.params=setTrapParams(list(barcode=c('AAA','BBB')))
#' ## not print
#' trap.params=setTrapParams(print=FALSE)
#' @export
setTrapParams <- function(trap.params=NULL, print=TRUE, check=TRUE, ...) {
  ## set default
  if (!is.null(trap.params)) {
    if (!is.list(trap.params)) stop("trap.params should be a list, see TRAP.PARAMS!")
    for (i in names(TRAP.PARAMS())) { # for pars not set in trap.params, use the default one.
      if (!(i %in% names(trap.params))) {
        trap.params[[i]]=TRAP.PARAMS()[[i]]
      }
    }
  } else {
    trap.params=TRAP.PARAMS() # default theme
  }

  ## delete non TRAP.PARAMS
  trap.params[!(names(trap.params) %in% names(TRAP.PARAMS()))]=NULL

  trap.params$tails.search=tolower(trap.params$tails.search)

  if (check) {
    if (!(trap.params$tails.search %in% c('genome','peaks','peak','no','none'))) stop("error tails.search, should be genome/peak/no")
  }
  if (trap.params$tails.search %in% c('peaks','peak')) trap.params$tails.search='peaks'
  if (trap.params$tails.search %in% c('no','none')) trap.params$tails.search='no'

  if (!print) return(trap.params)

  for (iname in names(TRAP.PARAMS())) {
    if (is.null(trap.params[[iname]])) {
      str=sprintf("%s = NULL", iname)
      .msg(string=str, showTime=FALSE, pre='', ...)
      next
    }
    i=which(names(trap.params)==iname)
    str=trap.params[[i]]
    if (is.character(trap.params[[i]])) {
      str=.shortStr(trap.params[[i]])
    }
    if (iname=='chrs') str=sprintf("%s (%d chrs)", str, length(trap.params[[i]]))
    if (iname=='barcode') str=sprintf("%s (%d barcodes)", str, length(trap.params[[i]]))
    str=sprintf("%s = %s", names(trap.params)[i], str)
    .msg(string=str, showTime=FALSE, pre='', ...)
  }

  return(trap.params)
}

#' Set tool names for running scAPAtrap
#'
#' setTools sets tool names with full path for running scAPAtrap and checks the validity.
#'
#' @param tools A list providing full path of four tools: samtools, umitools, featureCounts, star.
#' @param check Whether to check the validity of values (default is TRUE).
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{If `TRUE` basic status updates will be printed along the way.}
#' \item{logf }{If not NULL, then it should be a character string denoting a file name. Then message will be written to `logf`.}
#' }
#' @return A list with tool names as the names of the list.
#' @examples
#' tools=list(samtools='/home/dell/miniconda3/envs/scAPA/bin/samtools',
#'            umitools='/home/dell/miniconda3/envs/scAPA/bin/umi_tools',
#'            featureCounts="/home/dell/miniconda3/envs/scAPA/bin/featureCounts",
#'            star='/home/dell/miniconda3/envs/scAPA/bin/STAR')
#' \dontrun{
#' tools=setTools(tools)
#' }
#' tools=setTools(tools, check=FALSE)
#' @export
setTools <- function(tools, check=TRUE, ...) {
  TOOLS=c('samtools','umitools', 'featureCounts', 'star')

  tools[!(names(tools) %in% TOOLS)]=NULL

  if (check) {
    for (i in TOOLS) {
      if (!(i %in% names(tools))) stop(i, " not provided!")
    }
    for (i in tools) {
      if (!file.exists(i)) stop(i, " not exists!")
    }
  }

  for (i in 1:length(tools)) {
    str=sprintf("%s = %s", names(tools)[i], tools[[i]])
    .msg(string=str, showTime=FALSE, pre='', ...)
  }
  return(tools)
}

# .shortStr('hello')
# .shortStr(rep('chr', 200)) #--> [1] "chr,chr,chr,chr,chr,..."
.shortStr<-function(strs, len=20) {
  cs=paste0(strs, collapse=',')
  if (nchar(cs)<=len) return(cs)
  cs=paste0(substr(cs, 1, min(len, nchar(cs))), '...')
  return(cs)
}


## .printStep("step1")
## ...
## 2023-09-12 15:21:42 dedupByPos
## ...
.printStep <- function(str, ...) {
  .msg(string=paste0(rep('#', 80), collapse=''), showTime=FALSE, pre='', ...)
  .msg(string=str, showTime=TRUE, pre='', ...)
  .msg(string=paste0(rep('#', 80), collapse=''), showTime=FALSE, pre='', ...)
}

## filenames: if is NULL, then not check/stop (just print varname=NULL)
## varnames: if not NULL, then will print like "varname='filename'"
## .checkAndPrintFiles(c("/mnt/for.txt", "/mnt/rev.txt"), c("forwardPeaksFile", "reversePeaksFile")) ## will check and print
## .checkAndPrintFiles(NULL, c("forwardPeaksFile", "reversePeaksFile")) ## forwardPeaksFile=NULL\nreversePeaksFile=NULL
.checkAndPrintFiles <- function(filenames=NULL, varnames=NULL, ...) {

  nullFile=is.null(filenames)

  ## print output flies, similar to .printStep, but use '++++'
  if (!is.null(varnames)) {

    quote='\''
    if (length(filenames)!=length(varnames)) {
      filenames=rep('NULL', length(varnames))
      quote=''
    }

    .msg(string=paste0(rep('+', 80), collapse=''), showTime=FALSE, pre='', ...)

    for (i in 1:length(varnames)) {
      .msg(string=sprintf("%s=%s%s%s", varnames[i], quote, filenames[i], quote), showTime=FALSE, pre='', ...)
    }

    .msg(string=paste0(rep('+', 80), collapse=''), showTime=FALSE, pre='', ...)
  }

  ## check and stop
  if (!nullFile) {
    NO=FALSE
    for (f in filenames) {
      if (!file.exists(f)) {
        .msg(string=sprintf("%s not exists!", f), showTime=FALSE, pre='', ...)
        NO=TRUE
      }
    }
    if (NO) stop("scAPAtrap error, please check logf to see Logs (if set logf)!")
  }
  return(invisible(NULL))
}


## print TRAPFILES and note
.printTRAPFILES<-function(nc=50, note=TRUE, ...) {
  files=TRAPFILES()
  if (nrow(files)>0) nc=min(100, max(nchar(files$what)))
  if (note) .msg(string="Files marked with * can be savely deleted!", showTime=FALSE, pre='', ...)

  for (i in 1:nrow(files)) {
    str=sprintf(paste0("%", nc, "s ==> %s"), files[i, 2], files[i, 1])
    .msg(string=str, showTime=FALSE, pre='', ...)
  }
}




#' Initiate scAPAtrap and check/set parameters
#'
#' @inheritParams scAPAtrap
#' @return NULL. But the input tools, trap.params, inputBam, outputDir may be modifed during initiation.
#' @examples
#' \dontrun{
#' tools=list(samtools='/home/dell/miniconda3/envs/scAPA/bin/samtools',
#'            umitools='/home/dell/miniconda3/envs/scAPA/bin/umi_tools',
#'            featureCounts="/home/dell/miniconda3/envs/scAPA/bin/featureCounts",
#'            star='/home/dell/miniconda3/envs/scAPA/bin/STAR')
#' trap.params=TRAP.PARAMS()
#' trap.params$chr=1:3
#'
#' ## input BAM
#' dir0='/mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/xwu/dataFly/'
#' inputBam=paste0(dir0, "filterCB.bam")
#'
#' ## log file
#' logf=gsub('.bam', '.scAPAtrap.log', inputBam, fixed=TRUE)
#'
#' ## output dir (will be under inputBam's dir)
#' outputDir="APA.result"
#'
#' ## one step scAPAtrap
#' initScAPAtrap(tools=tools, trap.params=trap.params,
#'          inputBam=inputBam, outputDir=outputDir, logf=logf)
#' }
#' ## just test the function with toCheck=FALSE
#' initScAPAtrap(tools='just-for-fun.tools',
#'               trap.params=list('just-for-fun.params'),
#'               inputBam='just-for-fun.bam', outputDir='just-for-fun.dir', toCheck=FALSE)
#' trap.params
#' @export
initScAPAtrap<-function(tools, trap.params, inputBam, outputDir, ...)  {

  fn='initScAPAtrap'

  logf <- .getDotArg("logf", NULL, ...)
  if (!is.null(logf)) {
    if (file.exists(logf)) stop("logf is already exists, Please delete the LOG file or set another logf!")
    .msg(string=sprintf("LOG file = %s", logf), showTime=FALSE, pre='', ...)
  }

  ## toCheck=FALSE, just for debug test
  toCheck <- .getDotArg("toCheck", TRUE, ...)

  ## global=TRUE, to change global vars (used by step-by-step); otherwise return list of (trap.params... ) used by once-scAPAtrap
  global <- .getDotArg("global", TRUE, ...)

  .printStep('initiate scAPAtrap', ...)

  if (toCheck) .checkAndPrintFiles(inputBam, ...)

  if (!grepl('/', outputDir)) {
    outputDir = paste0(dirname(inputBam), '/', outputDir)
  }

  if (toCheck) {
    if (file.exists(outputDir)) stop('initScAPAtrap error: outputDir [', outputDir, '] exists, please change another dirname to avoid overwriting\n')
  }

  .msg(string=sprintf("inputBam = %s", inputBam), showTime=FALSE, pre='', ...)
  .msg(string=sprintf("outputDir = %s", outputDir), showTime=FALSE, pre='', ...)

  trap.params <- setTrapParams(trap.params, print=FALSE, check = toCheck, ...)
  tools <- setTools(tools, check = toCheck, ...)

  TRAPFILES(clear=toCheck)

  ## get chr names: if not provided, then get from BAM files
  if (!toCheck) {
    # just print
    x=setTrapParams(trap.params, print=TRUE, check = FALSE, ...)
    return(invisible(NULL))
  }

  if (length(trap.params$chrs)==0) {
      trap.params$chrs <- .getBAMchrs(inputBam)
    .msgin(fn='', string=sprintf("Total %d chrs (e.g., %s) from inputBam", length(trap.params$chrs), .shortStr(trap.params$chrs, 50)), ...)

    if (length(trap.params$chrs)>100) {
      .msgin(fn='', string=sprintf("Too many chrs (%d, >100) from inputBam, please use trap.params$chrs to explicitly set chrs", length(trap.params$chrs)), ...)
      stop("initScAPAtrap quit because of too many chrs from inputBam, please check logf to see Logs!")
    }
  }
  if (length(trap.params$chrs)==0) {
    .msg(string="No chrs in inputBam or in trap.params!", showTime=FALSE, pre='', ...)
    stop("initScAPAtrap error because of no chrs provided, please check logf to see Logs!")
  }

  ## just print
  trap.params <- setTrapParams(trap.params, print=TRUE, check = toCheck, ...)

  if (global) {
    trap.params <<- trap.params
    tools <<- tools
    inputBam <<- inputBam
    outputDir <<- outputDir
    return(invisible(NULL))
  } else {
    return(list(trap.params=trap.params, tools=tools, inputBam=inputBam, outputDir=outputDir))
  }

}

#' Wrapper function for one-step running of scAPAtrap
#'
#' @param tools A list contains four tools used in scAPAtrap.
#' @param trap.params Parameters of running scAPAtrap, see global variable \code{\link{TRAP.PARAMS}}.
#' @param inputBam Input bam file name. If trap.params$chr is NULL, will try to get chrs from inputBam.
#' @param outputDir Output dir for storing scAPAtrap's final output -- peaks and counts files.
#' If only dirname is given (e.g., APAres but not ./APAres), then will add path of inputBam to the dir name.
#' If outputDir already exists, an error will be raised to avoid overwritting.
#' However, many temporary files are generated during running scAPAtrap, which could be seen from global variable TRAPFILES. Those files marked with * could be savely deleted.
#' @param ... Arguments passed to other methods and/or advanced arguments.
#' Advanced arguments:
#' \describe{
#' \item{verbose }{If `TRUE` basic status updates will be printed along the way.}
#' \item{logf }{If not NULL, then it should be a character string denoting a file name.
#' If logf is a character string, then message will be written to `logf`.
#' But if logf already exists, this function will quit to avoid overwriting an existing log file.
#' The logf will log full information during scAPAtrap, including time, command, and output files.}
#' }
#' @return A file name storing scAPAtrapData object, which is <outputDir>/scAPAtrapData.rda
#' @examples
#' \dontrun{
#' tools=list(samtools='/home/dell/miniconda3/envs/scAPA/bin/samtools',
#'            umitools='/home/dell/miniconda3/envs/scAPA/bin/umi_tools',
#'            featureCounts="/home/dell/miniconda3/envs/scAPA/bin/featureCounts",
#'            star='/home/dell/miniconda3/envs/scAPA/bin/STAR')
#' trap.params=TRAP.PARAMS()
#' trap.params$chrs=c('2L','2R','3L','3R','4','X','Y')
#'
#' ## input BAM
#' dir0='/mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/xwu/dataFly/'
#' inputBam=paste0(dir0, "filterCB.bam")
#'
#' ## log file
#' logf=gsub('.bam', '.scAPAtrap.log', inputBam, fixed=TRUE)
#' unlink(logf)
#'
#' ## barcode (if have)
#' barcodefile <- paste0(dir0, 'barcode.17.txt')
#' barcode <- utils::read.delim2(barcodefile, header = FALSE)
#' barcode <- gsub('-[0-9]','',barcode$V1)
#' ## add barcode to trap.params
#' trap.params$barcode=barcode
#'
#' ## output dir (will be under inputBam's dir)
#' outputDir="APA.result"
#'
#' ## one step scAPAtrap
#' scAPAtrap(tools=tools, trap.params=trap.params,
#'          inputBam=inputBam, outputDir=outputDir, logf=logf)
#' }
#' @export
scAPAtrap <- function(tools, trap.params, inputBam, outputDir, ...) {

  #set.seed(123)

  ################### 1. setting #######################
  ## set logf, trap.params, tools, TRAPFILES, chrs (if is NULL)
  init=initScAPAtrap(tools=tools, trap.params=trap.params, inputBam=inputBam, outputDir=outputDir, global=FALSE, ...)
  tools=init$tools
  trap.params=init$trap.params
  inputBam=init$inputBam
  outputDir=init$outputDir

  ################## 2. findUniqueMap #######################
  # findUniqueMap: uniq >> sort >> index with samtools
  # Skip this step if the inputBam contains only unique mappings.
  if (trap.params$findUniqueMap) {
    .printStep('findUniqueMap', ...)
    inputBam.u <- findUniqueMap(tools$samtools, inputBam, thread=trap.params$thread,
                                sort = trap.params$findUniqueMap, index = trap.params$findUniqueMap.index, ...)
  } else {
    inputBam.u=inputBam
  }
  .checkAndPrintFiles(filenames = inputBam.u, varnames=c("inputBam.u"), ...)

  #################### 3. dedupByPos #####################
  ## dedupByPos: remove duplicates with umitools
  ## If the inputBAM is from non-10X protocols, set TenX=FALSE
  ## However, if STARsolo is used for CEL-seq which contains cell barcode tag (CB) and UMI barcode tag (UB), set TenX=FALSE
  .printStep('dedupByPos', ...)
  inputBam.u.dedup <- dedupByPos(tools$umitools, input=inputBam.u, TenX=trap.params$TenX, ...)
  .checkAndPrintFiles(filenames = inputBam.u.dedup, varnames=c("inputBam.u.dedup"), ...)

  #################### 4. separateBamBystrand #####################
  ## separateBamBystrand: split the BAM file by strand to forward.bam and reverse.bam
  .printStep('separateBamBystrand', ...)
  inputBam.u.dedup.seps <- separateBamBystrand(tools$samtools, input=inputBam.u.dedup, thread=trap.params$thread, ...)
  .checkAndPrintFiles(filenames = inputBam.u.dedup.seps, varnames=c("inputBam.u.dedup.seps[1]", "inputBam.u.dedup.seps[2]"), ...)

  ################### 5. findPeaksByStrand ######################
  ## findPeaksByStrand: find peaks based on BAM coverages
  ## This step merges previous loadBpCoverages and findPeaks.
  .printStep('findPeaksByStrand', ...)
  forwardPeaksFile=paste0(inputBam.u.dedup.seps[1], '.peaks')
  forwardPeaksFile <-findPeaksByStrand(bamFile=inputBam.u.dedup.seps[1], chrs=trap.params$chrs, strand='+',
                                   L=trap.params$readlength, maxwidth=trap.params$maxwidth,
                                   cutoff = trap.params$cov.cutoff,
                                   ofile=forwardPeaksFile, ...)
  .checkAndPrintFiles(filenames = forwardPeaksFile, varnames=c("forwardPeaksFile"), ...)

  reversePeaksFile=paste0(inputBam.u.dedup.seps[2], '.peaks')
  reversePeaksFile <-findPeaksByStrand(bamFile=inputBam.u.dedup.seps[2], chrs=trap.params$chrs, strand='-',
                                   L=trap.params$readlength, maxwidth=trap.params$maxwidth,
                                   cutoff = trap.params$cov.cutoff,
                                   ofile=reversePeaksFile, ...)
  .checkAndPrintFiles(filenames = reversePeaksFile, varnames=c("reversePeaksFile"), ...)

  ################## 6. generateSAF #######################
  ## generateSAF: generate a .SAF file that records the information of identified peaks
  .printStep('generateSAF', ...)
  peaksfile <- generateSAF(forwardPeaks=forwardPeaksFile, reversePeaks=reversePeaksFile, outputdir=outputDir, ...)
  .checkAndPrintFiles(filenames = peaksfile, varnames=c("peaksfile"), ...)

  ################# 7. generateFinalBam ########################
  ## generateFinalBam: get a small BAM that contains the peak regions.
  ## Here the inputBam is used, which will contain duplicated reads.
  .printStep('generateFinalBam', ...)
  final.bam <- generateFinalBam(tools$featureCounts, tools$samtools, input=inputBam, peakfile=peaksfile, thread=trap.params$thread, ...)
  .checkAndPrintFiles(filenames = final.bam, varnames=c('final.bam'), ...)

  ################# 8. countPeaks ########################
  ## countPeaks: count read number in each peak with umitools
  .printStep('countPeaks', ...)
  countsfile<- countPeaks(tools$umitools, input=final.bam, outputdir=outputDir, TenX=trap.params$TenX, ...)
  .checkAndPrintFiles(filenames = countsfile, varnames=c('countsfile'), ...)
  countsPeaksFiles=list(countsfile=countsfile, peaksfile=peaksfile)

  ################## 9. findTails #######################
  ## findTails: peaks, genome, no
  if (trap.params$tails.search=='peaks') {
    .printStep('reducePeaks (when tails.search=peaks)', ...)
    #8.1# reducePeaks: remove small peaks to speed up tail-searching
    countsPeaksFiles <- reducePeaks(countsfile=countsfile, peaksfile=peaksfile,
                                         min.cells = trap.params$min.cells, min.count = trap.params$min.count,
                                         suffix='.reduced', toSparse=FALSE, ...)
    peaksfile.reduced=countsPeaksFiles$peaksfile
    .checkAndPrintFiles(filenames = peaksfile.reduced, varnames=c('peaksfile.reduced'), ...)

    .printStep('findTailsByPeaks (when tails.search=peaks)', ...)
    #8.2# findTailsByPeaks: find tails nearing peaks. Here the inputBam is used, which will contain duplicated reads.
    tailsfile=paste0(inputBam, '.peaks.tails')
    tailsfile <- findTailsByPeaks(bamfile = inputBam, peaksfile=peaksfile.reduced, d=trap.params$d, tailsfile=tailsfile, ...)

  } else if (trap.params$tails.search=='genome') {
    .printStep('findTails (when tails.search=genome)', ...)
    #8# findTails: find tails genome-wide. Here the inputBam is used, which will contain duplicated reads.
    tailsfile=paste0(inputBam, '.genome.tails')
    tailsfile <- findTails(bamfile = inputBam, tailsfile=tailsfile)

  } else {
    ## not find tails
    tailsfile=NULL
  }
  .checkAndPrintFiles(filenames = tailsfile, varnames=c('tailsfile'), ...)
  ################### 10. generatescExpMa ######################
  ## generatescExpMa: generate final PA data including peaks and counts information
  .printStep('generatescExpMa', ...)
  outputfile=paste0(outputDir, '/scAPAtrapData.rda')
  outputfile <- generatescExpMa(countsfile=countsPeaksFiles$countsfile,
                            peaksfile=countsPeaksFiles$peaksfile,
                            barcode=trap.params$barcode,
                            tails=tailsfile, d=trap.params$d,
                            min.cells=trap.params$min.cells, min.count=trap.params$min.count,
                            ofile=outputfile, ...)
  .checkAndPrintFiles(filenames = outputfile, varnames=c('outputfile'), ...)

  ################### 11. clean TRAPFILES ######################
  ## show TRAPFILES for cleanning
  .printStep('All output files in running order', ...)
  .printTRAPFILES(nc=50, ...)

  ## return scAPAtrapData.rda
  return(outputfile)

}





