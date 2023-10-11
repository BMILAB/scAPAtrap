## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 6,
  fig.height = 5.5,
  collapse = TRUE,
  warning = FALSE,
  comment = "#>"
)

## ----install, eval=FALSE------------------------------------------------------
#  devtools::install_github("BMILAB/scAPAtrap")

## ----one_step_scAPAtrap, eval=FALSE-------------------------------------------
#  library(scAPAtrap)
#  
#  ## input BAM
#  dir0='/demoFly/'
#  inputBam=paste0(dir0, "fly_demo.bam")
#  
#  ## log file to LOG all information (time, command, output file names..)
#  logf=gsub('.bam', '.APA.notails.onestep.log', inputBam, fixed=TRUE)
#  
#  ## output dir (will be under inputBam's dir if only dirname is provided; otherwise use the full path)
#  outputDir="APA.notails.onestep"
#  
#  ## full path of tools
#  tools=list(samtools='/home/dell/miniconda3/envs/scAPA/bin/samtools',
#             umitools='/home/dell/miniconda3/envs/scAPA/bin/umi_tools',
#             featureCounts="/home/dell/miniconda3/envs/scAPA/bin/featureCounts",
#             star='/home/dell/miniconda3/envs/scAPA/bin/STAR')
#  
#  ## set parameters, first load default parameters and then modify
#  trap.params=setTrapParams()
#  ## set tail search way
#  trap.params$tails.search='no'
#  ## set chromosome names
#  trap.params$chrs=c('2L','2R','3L','3R','4','X','Y')
#  
#  ## we can also get chromosome names from the BAM file
#  ## however, we should also check the number of chrs to avoid too many sccafolds or fragments.
#  chrs=scAPAtrap:::.getBAMchrs(inputBam)
#  length(chrs)
#  
#  ## barcode (if have)
#  trap.params$barcode <- read.delim2(paste0(dir0, 'barcode.txt'), header = F)$V1
#  
#  ## Run scAPAtrap
#  scAPAtrap(tools=tools,
#            trap.params=trap.params,
#            inputBam=inputBam,
#            outputDir=outputDir,
#            logf=logf)
#  

## ----step1_initiation, eval=FALSE---------------------------------------------
#  ################### 1. initiate #######################
#  library(scAPAtrap)
#  
#  ## input BAM and output dir
#  dir0='/mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/xwu/demoFly/'
#  inputBam=paste0(dir0, "fly_demo.bam")
#  outputDir='APA.notails.stepbystep'
#  
#  ## logf
#  logf=gsub('.bam', '.APA.notails.stepbystep.log', inputBam, fixed=TRUE)
#  unlink(logf)
#  
#  ## tools
#  tools=list(samtools='/home/dell/miniconda3/envs/scAPA/bin/samtools',
#             umitools='/home/dell/miniconda3/envs/scAPA/bin/umi_tools',
#             featureCounts="/home/dell/miniconda3/envs/scAPA/bin/featureCounts",
#             star='/home/dell/miniconda3/envs/scAPA/bin/STAR')
#  
#  ## trap.params: first get default parameters and set specific ones
#  trap.params=setTrapParams()
#  
#  ## chrs
#  trap.params$chrs=c('2L','2R','3L','3R','4','X','Y')
#  
#  ## barcode (if have)
#  trap.params$barcode <- read.delim2(paste0(dir0, 'barcode.txt'), header = F)$V1
#  
#  ## initiation, will print these settings in logf and on screen
#  initScAPAtrap(tools=tools,
#                trap.params=trap.params,
#                inputBam=inputBam,
#                outputDir=outputDir,
#                logf=logf)

## ----step2_findUniqueMap, eval=FALSE------------------------------------------
#  ################## 2. findUniqueMap #######################
#  # findUniqueMap: uniq >> sort >> index with samtools
#  # Skip this step if the inputBam contains only unique mappings.
#  inputBam.u=inputBam
#  if (trap.params$findUniqueMap)
#    inputBam.u <- findUniqueMap(tools$samtools,
#                                inputBam,
#                                thread=trap.params$thread,
#                                sort = trap.params$findUniqueMap,
#                                index = trap.params$findUniqueMap.index,
#                                logf=logf)
#  
#  scAPAtrap:::.checkAndPrintFiles(varnames=c("inputBam.u"),
#                                  filenames=inputBam.u, logf=logf)

## ----step3_dedupByPos, eval=FALSE---------------------------------------------
#  #################### 3. dedupByPos #####################
#  ## dedupByPos: remove duplicates with umitools
#  inputBam.u.dedup <- dedupByPos(tools$umitools,
#                                 input=inputBam.u,
#                                 TenX=trap.params$TenX,
#                                 logf=logf)
#  
#  scAPAtrap:::.checkAndPrintFiles(varnames=c("inputBam.u.dedup"),
#                                  filenames=inputBam.u.dedup, logf=logf)

## ----step4_separateBamBystrand, eval=FALSE------------------------------------
#  #################### 4. separateBamBystrand #####################
#  ## separateBamBystrand: split the BAM file by strand to forward.bam and reverse.bam
#  inputBam.u.dedup.seps <- separateBamBystrand(tools$samtools,
#                                               input=inputBam.u.dedup,
#                                               thread=trap.params$thread,
#                                               logf=logf)
#  
#  scAPAtrap:::.checkAndPrintFiles(varnames=c("inputBam.u.dedup.seps[1]", "inputBam.u.dedup.seps[2]"),
#                                  filenames=inputBam.u.dedup.seps)

## ----step5_findPeaksByStrand, eval=FALSE--------------------------------------
#  ################### 5. findPeaksByStrand ######################
#  ## findPeaksByStrand: find peaks based on BAM coverages
#  forwardPeaksFile=paste0(inputBam.u.dedup.seps[1], '.peaks')
#  forwardPeaksFile <-findPeaksByStrand(bamFile=inputBam.u.dedup.seps[1],
#                                       chrs=trap.params$chrs,
#                                       strand='+',
#                                       L=trap.params$readlength,
#                                       maxwidth=trap.params$maxwidth,
#                                       cutoff = trap.params$cov.cutoff,
#                                       ofile=forwardPeaksFile,
#                                       logf=logf)
#  
#  scAPAtrap:::.checkAndPrintFiles(varnames='forwardPeaksFile',
#                                  filenames=forwardPeaksFile, logf=logf)
#  
#  reversePeaksFile=paste0(inputBam.u.dedup.seps[2], '.peaks')
#  reversePeaksFile <-findPeaksByStrand(bamFile=inputBam.u.dedup.seps[2],
#                                       chrs=trap.params$chrs, strand='-',
#                                       L=trap.params$readlength,
#                                       maxwidth=trap.params$maxwidth,
#                                       cutoff = trap.params$cov.cutoff,
#                                       ofile=reversePeaksFile, logf=logf)
#  
#  scAPAtrap:::.checkAndPrintFiles(varnames='reversePeaksFile',
#                                  filenames=reversePeaksFile, logf=logf)

## ----step6_generateSAF, eval=FALSE--------------------------------------------
#  ################## 6. generateSAF #######################
#  ## generateSAF: generate a .SAF file that records the information of identified peaks
#  peaksfile <- generateSAF(forwardPeaks=forwardPeaksFile,
#                           reversePeaks=reversePeaksFile,
#                           outputdir=outputDir,
#                           logf=logf)
#  
#  scAPAtrap:::.checkAndPrintFiles(varnames='peaksfile',
#                                  filenames=peaksfile, logf=logf)

## ----step7_generateFinalBam, eval=FALSE---------------------------------------
#  ################# 7. generateFinalBam ########################
#  ## generateFinalBam: get a small BAM that contains the peak regions.
#  ## Here the inputBam is used, which will contain duplicated reads.
#  final.bam <- generateFinalBam(tools$featureCounts,
#                                tools$samtools,
#                                input=inputBam,
#                                peakfile=peaksfile,
#                                thread=trap.params$thread,
#                                logf=logf)
#  
#  scAPAtrap:::.checkAndPrintFiles(varnames='final.bam',
#                                  filenames=final.bam, logf=logf)

## ----step8_countPeaks, eval=FALSE---------------------------------------------
#  ################# 8. countPeaks ########################
#  ## countPeaks: count read number in each peak with umitools
#  countsfile<- countPeaks(tools$umitools,
#                          input=final.bam,
#                          outputdir=outputDir,
#                          TenX=trap.params$TenX,
#                          logf=logf)
#  
#  scAPAtrap:::.checkAndPrintFiles(varnames='countsfile',
#                                  filenames=countsfile, logf=logf)
#  
#  ## record output file for reducePeaks (if tails.search=peaks)
#  

## ----step9_findTails_no, eval=FALSE-------------------------------------------
#  ################## 9. findTails #######################
#  ## tails.search=no
#  tailsfile=NULL
#  
#  countsPeaksFiles=list(countsfile=countsfile, peaksfile=peaksfile)
#  
#  scAPAtrap:::.checkAndPrintFiles(varnames='tailsfile',
#                                  filenames=tailsfile, logf=logf)

## ----step9.2_findTails_genome_peaks, eval=FALSE-------------------------------
#  
#  ## findTails: peaks, genome
#  
#  ## tails.search=peaks
#  if (trap.params$tails.search=='peaks') {
#    #8.1# reducePeaks: remove small peaks to speed up tail-searching
#    countsPeaksFiles <- reducePeaks(countsfile=countsfile,
#                                    peaksfile=peaksfile,
#                                    min.cells = trap.params$min.cells,
#                                    min.count = trap.params$min.count,
#                                    suffix='.reduced', logf=logf)
#  
#    peaksfile.reduced=countsPeaksFiles$peaksfile
#  
#    scAPAtrap:::.checkAndPrintFiles(varnames='peaksfile.reduced',
#                                    filenames=peaksfile.reduced, logf=logf)
#  
#    #8.2# findTailsByPeaks: find tails nearing peaks.
#    tailsfile=paste0(inputBam, '.peaks.tails')
#    tailsfile <- findTailsByPeaks(bamfile = inputBam,
#                                  peaksfile=peaksfile.reduced,
#                                  d=trap.params$d,
#                                  tailsfile=tailsfile, logf=logf)
#  
#    scAPAtrap:::.checkAndPrintFiles(varnames='tailsfile',
#                                    filenames=tailsfile, logf=logf)
#  
#  ## tails.search=genome
#  } else if (trap.params$tails.search=='genome') {
#    #8# findTails: find tails genome-wide.
#    tailsfile=paste0(inputBam, '.genome.tails')
#    tailsfile <- findTails(bamfile = inputBam,
#                           tailsfile=tailsfile, logf=logf)
#  
#    scAPAtrap:::.checkAndPrintFiles(varnames='tailsfile',
#                                    filenames=tailsfile, logf=logf)
#  }
#  

## ----step10_generatescExpMa, eval=FALSE---------------------------------------
#  ################### 10. generatescExpMa ######################
#  ## generatescExpMa: generate final PA data including peaks and counts information
#  outputfile=paste0(outputDir, '/scAPAtrapData.rda')
#  outputfile <- generatescExpMa(countsfile=countsPeaksFiles$countsfile,
#                                peaksfile=countsPeaksFiles$peaksfile,
#                                barcode=trap.params$barcode,
#                                tails=tailsfile,
#                                d=trap.params$d,
#                                min.cells=trap.params$min.cells,
#                                min.count=trap.params$min.count,
#                                ofile=outputfile,
#                                logf=logf)
#  ## Final scAPAtrapData.rda
#  scAPAtrap:::.checkAndPrintFiles(varnames='outputfile',
#                                  filenames=outputfile, logf=logf)

## ----step11_clean, eval=FALSE-------------------------------------------------
#  ################### 10. clean TRAPFILES ######################
#  ## show TRAPFILES for cleanning
#  scAPAtrap:::.printTRAPFILES(nc=50, logf=logf)

## ----movAPA, eval=FALSE-------------------------------------------------------
#  install.packages("devtools")
#  require(devtools)
#  install_github("BMILAB/movAPA")
#  library(movAPA)
#  browseVignettes('movAPA')

## ----vizAPA, eval=FALSE-------------------------------------------------------
#  install_github("BMILAB/vizAPA")
#  library(vizAPA)
#  browseVignettes('vizAPA')

## ----search_tails_separately, eval=FALSE--------------------------------------
#  
#  #### to search tails genome-wide
#  tailsfile=paste0(inputBam, '.genome.tails')
#  tailsfile <- findTails(bamfile = inputBam, tailsfile=tailsfile)
#  
#  #### or to search tails peak-wide
#  # After running the scAPAtrap's pipeline, countsfile and peaksfile stores the peaks.meta and peak.counts.
#  # It is easy to filter peaks that are not highly-expressed but moderately expressed.
#  
#  # first, generate a peak file with reduced peaks, here peaks with moderately expression level can be used
#  # For example, if we consider peaks with 50 cells and 100 counts are definitely true peaks,
#  # and consider peaks with <10 cells and <20 counts are definitely false peaks,
#  # but peaks between 10~50 cells and 20 can 100 counts could be true or false.
#  # Then we can retain these moderately peaks for tail searching.
#  
#  # to retain peaks with >=10 cells but <50 cells and >=20 counts and <100 cells as unsure peaks
#  peaks.unsure <- reducePeaks(countsfile, peaksfile,
#                          min.cells = 10, min.count = 20,
#                          max.cells = 49, max.count = 99,
#                          suffix='.reduced', ...)
#  
#  # to retain those definitely true peaks too (>=50 cells, and >=100 counts)
#  peaks.true <- reducePeaks(countsfile, peaksfile,
#                          min.cells = 50, min.count = 100,
#                          max.cells = NULL, max.count = NULL,
#                          suffix='.true', ...)
#  
#  # to search tails for peaks.unsure
#  tailsfile=paste0(inputBam, '.peaks.unsure.tails')
#  tailsfile <- findTailsByPeaks(bamfile = inputBam, peaksfile=peaks.unsure, d=trap.params$d, tailsfile=tailsfile, logf=logf)
#  
#  # then use these tails to filter peaks, which are unsure-peaks supported by tails
#  outputfile=paste0(tailsfile,'.peaks.rda')
#  outputfile <- generatescExpMa(countsfile=peaks.unsure[1],
#                                peaksfile=peaks.unsure[2],
#                                barcode=trap.params$barcode,
#                                tails=tailsfile, d=trap.params$d,
#                                min.cells=1, min.count=1,
#                                ofile=outputfile, logf=logf)
#  
#  # then we can combine peaks.true and outputfile
#  
#  

## ----change_chr_names, eval=FALSE---------------------------------------------
#  ## some code like this to remove 'chr'
#  dataframe$chrs=gsub('chr', '', dataframe$chrs)
#  
#  ## some code like this to add 'chr'
#  dataframe$chrs=paste0('chr', dataframe$chrs)

