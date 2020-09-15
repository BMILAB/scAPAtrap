# 2020/7/18 funclibs for scAPAtrap
#-----------------------------------Preprocess----------------------------------------------
#' Filter the BAM file to preserve the unique comparison results
#'
#' @param samtools.path The path of the samtools
#' @param input BAM file
#' @param thread Number of threads
#' @param sort Logical value, whether to sort
#' @param index Logical value, whether or not to build the index
#' @examples
#' samtools.path <- '/home/aa/miniconda2/envs/nar_env/bin/samtools'
#' input <- './data/demo.bam'
#' findUniqueMap(samtools.path,input,24)
#' @return Output path of the calculated BAM file
#' @export
findUniqueMap <- function(samtools.path, input, thread, sort = T, index = T) {
    output <- paste0(gsub(".bam$", "", input), ".Uniq.bam")
    command <- paste0(samtools.path, " view -@ ", thread, " -h -F 256 -bS ", input, " > ", output)
    system(command = command, wait = T)

    res <- output
    if (sort) {
        output.sorted <- paste0(gsub(".Uniq.bam$", "", output), ".UniqSorted.bam")
        command <- paste0(samtools.path, " sort -@ ", thread, " -o ", output.sorted, " ", output)
        system(command = command, wait = T)
        res <- output.sorted
        file.remove(output)
    }
    if (index) {
        command <- paste0(samtools.path, " index -@ ", thread, " ", res)
        system(command = command, wait = T)
    }

    return(res)
}


#' Deduplication according to the comparison position
#'
#' @param umitools.path The path of the umi_tools
#' @param input BAM file,but need index files
#' @param TenX Logical value. Is it 10x data?
#' @examples
#' umitools.path <- '/home/aa/miniconda2/envs/umi_tools/bin/umi_tools'
#' input <- './data/demo.Uniq.sorted.bam'
#' dedupByPos(umitools.path,input)
#' @return Output path of the calculated BAM file
#' @export
dedupByPos <- function(umitools.path, input, TenX = T) {
    output <- paste0(gsub(".bam$", "", input), ".dedup.bam")
    if (TenX) {
        dedup.command <- paste0(umitools.path, " dedup -I ", input, " -S ", output, " --method=unique --extract-umi-method=tag --umi-tag=UB --cell-tag=CB")
        system(command = dedup.command, wait = T)
    } else {
        dedup.command <- paste0(umitools.path, " dedup -I ", input, " -S ", output, " --method=unique")
        system(command = dedup.command, wait = T)
    }
    file.remove(input)
    file.remove(paste0(input, ".bai"))
    return(output)
}


#' Separate BAM files according to strand
#'
#'@param samtools.path The path of the samtools
#'@param input BAM file
#'@param thread Number of threads
#'@examples
#'samtools.path <- '/home/aa/miniconda2/envs/nar_env/bin/samtools'
#'input <- './data/demo.Uniq.sorted.dedup.bam'
#'separateBamBystrand(samtools.path,input,24)
#'@return String vector, the path of the forward BAM file and the reverse BAM file generated after separating the BAM
#'@export
separateBamBystrand <- function(samtools.path, input, thread) {
    outputF <- paste0(gsub(".bam", "", input), ".forward.bam")
    outputR <- paste0(gsub(".bam", "", input), ".reverse.bam")
    samtools.Fcommand <- paste0(samtools.path, " view -@ ", thread, " -h -F 0x10 -bS ", input, " > ", outputF)
    system(command = samtools.Fcommand, wait = T)
    samtools.Rcommand <- paste0(samtools.path, " view -@ ", thread, " -h -f 0x10 -bS ", input, " > ", outputR)
    system(command = samtools.Rcommand, wait = T)

    command <- paste0(samtools.path, " index -@ ", thread, " ", outputF)
    system(command = command, wait = T)
    command <- paste0(samtools.path, " index -@ ", thread, " ", outputR)
    system(command = command, wait = T)

    file.remove(input)
    return(c(outputF, outputR))
}


#' Extract barcode and UMI in sequence
#'
#' @param umitools.path The path of the umi_tools
#' @param pattern Mode used to represent barcode and UMI (e.g The barcode length is 8 and the UMI length is 4, then pattern=CCCCCCCCNNNN)
#' @param input.seq path of fastq file,must contain two paths, the fastq of the first path records barcode and umi information
#' @param whitelist whitelist barcode
#' @examples
#' input.seq <- c('./fastq/SRR1610598_1.fastq','./fastq/SRR1610598_2.fastq')
#' pattern <- 'CCCCCCCCNNNN'
#' whitelist <- './fastq/whitelist.txt'
#' extractBcAndUb(umitools.path,pattern,input.seq,whitelist)
#' @export
extractBcAndUb <- function(umitools.path, pattern, input.seq, whitelist) {
    output1 <- paste0(gsub(".fastq$", "", input.seq[1]), ".extracted.fastq")
    output2 <- paste0(gsub(".fastq$", "", input.seq[2]), ".extracted.fastq")
    command <- paste0(umitools.path, " extract --bc-pattern=", pattern, " --stdin ", input.seq[1], " --stdout ", output1,
                      " --read2-in ", input.seq[2], " --read2-out=", output2, " --filter-cell-barcode ", " --whitelist=", whitelist)

    system(command = command, wait = T)

    return(output2)
}


#' Construct reference genome index
#'
#' @param star.path The path of the star.path
#' @param genome.fasta Reference genome sequence, fasta format
#' @param indexdir Directory of Reference Genome Index
#' @param thread Number of threads
#' @return Directory of Reference Genome Index
#' @examples
#' genome.fasta <- '/home/aa/genome_data/mm10/Mus_musculus.GRCm38.dna.primary_assembly.fa'
#' indexdir <- './mm10_index/'
#' generateRefIndex(star.path,genome.fasta,indexdir,24)
#' @export
generateRefIndex <- function(star.path, genome.fasta, indexdir, thread) {
    command <- paste0(star.path, " --runMode genomeGenerate --genomeDir ", indexdir, " --genomeFastaFiles ", genome.fasta,
                      " --runThreadN ", thread)

    system(command = command, wait = T)

    return(indexdir)
}



#' Align the sequence to the reference genome
#'
#' @param star.path The path of the star.path
#' @param indexdir Directory of Reference Genome Index
#' @param prefix The prefix of the output file name
#' @param thread Number of threads
#' @examples
#' input.seq <- 'SRR1610598_2.extracted.fastq'
#' generateAlignBam(star.path,indexdir,input.seq,'./BAM/rep1',12)
#' @return The path of the bam file generated after alignment
#' @export
generateAlignBam <- function(star.path, indexdir, input.seq, prefix, thread) {
    command <- paste0(star.path, " --runThreadN ", thread, " --genomeDir ", indexdir, " --readFilesIn ", input.seq, " --outFilterMultimapNmax 1 ",
                      " --outSAMtype BAM SortedByCoordinate ", " --outFileNamePrefix ", prefix)
    system(command = command, wait = T)
}

#-----------------------------------findTails----------------------------------------------
.readChrInfo <- function(file) {
    chr_length <- read.delim(file = file)
    colnames(chr_length) <- c("flag", "chr", "coord")
    chr_length$chr <- gsub("SN:", "", chr_length$chr)
    chr_length$coord <- as.numeric(gsub("LN:", "", chr_length$coord))
    chr_length <- subset(chr_length, chr != "MT")
    chr.g <- with(chr_length, GRanges(seqnames = chr, ranges = IRanges(start = 1, end = coord)))

    return(chr.g)
}

#' Find the precise polyA site on each chromosome
#'
#' @description Traverse each chromosome and find the site through the A-rich at the end. It should be noted that BAM files need to be indexed
#' @param bamfile BAM file
#' @param chr chromosome
#' @param len Chromosome length
#' @examples
#' findChrTails(bamfile,1,1000000001)
#' @export
findChrTails <- function(bamfile, chr, len) {
    which <- GenomicRanges::GRanges(seqnames = chr,start=1,end=len)
    what <- c("rname", "strand", "pos", "cigar", "seq")
    param <- Rsamtools::ScanBamParam(what = what, which = which)
    gal1 <- GenomicAlignments::readGAlignments(bamfile, use.names = TRUE, param = param)
    s_1 <- (grepl("[0-9]*M[1-9]{2,}S", gal1@cigar) & as.vector(gal1@strand) == "+")
    s_2 <- (grepl("[1-9]{2,}S[0-9]*M", gal1@cigar) & as.vector(gal1@strand) == "-")

    bam1 <- gal1[s_1]
    bam2 <- gal1[s_2]

    bam1 <- bam1[grepl("(A{3,}[^A]{0,2})*A{6,}([^A]{0,2}A{3,})*.{0,2}?$", bam1@elementMetadata@listData$seq)]
    bam2 <- bam2[grepl("^.{0,2}?(T{3,}[^T]{0,2})*T{6,}([^T]{0,2}T{3,})*", bam2@elementMetadata@listData$seq)]

    final_bam1 <- data.frame(chr = as.vector(seqnames(bam1)), strand = as.vector(strand(bam1)), coord = end(bam1))
    final_bam2 <- data.frame(chr = as.vector(seqnames(bam2)), strand = as.vector(strand(bam2)), coord = start(bam2))

    bam <- rbind(final_bam1, final_bam2)
    bam <- dplyr::group_by(bam, chr, strand, coord) %>% summarise(count = n())

    return(bam)
}


#' Find the precise polyA site
#'
#' @param bamfile BAM file
#' @export
findTails <- function(bamfile) {
    what <- c("rname", "strand", "pos", "cigar", "seq")
    param <- Rsamtools::ScanBamParam(what = what)
    gal1 <- GenomicAlignments::readGAlignments(bamfile, use.names = TRUE, param = param)
    s_1 <- (grepl("[0-9]*M[1-9]{2,}S", gal1@cigar) & as.vector(gal1@strand) == "+")
    s_2 <- (grepl("[1-9]{2,}S[0-9]*M", gal1@cigar) & as.vector(gal1@strand) == "-")

    bam1 <- gal1[s_1]
    bam2 <- gal1[s_2]

    bam1 <- bam1[grepl("(A{3,}[^A]{0,2})*A{6,}([^A]{0,2}A{3,})*.{0,2}?$", bam1@elementMetadata@listData$seq)]
    bam2 <- bam2[grepl("^.{0,2}?(T{3,}[^T]{0,2})*T{6,}([^T]{0,2}T{3,})*", bam2@elementMetadata@listData$seq)]

    final_bam1 <- data.frame(chr = as.vector(seqnames(bam1)), strand = as.vector(strand(bam1)), coord = end(bam1))
    final_bam2 <- data.frame(chr = as.vector(seqnames(bam2)), strand = as.vector(strand(bam2)), coord = start(bam2))

    bam <- rbind(final_bam1, final_bam2)
    bam <- dplyr::group_by(bam, chr, strand, coord) %>% summarise(count = n())

    chr <- colnames(bam)[1]
    return(bam)
}

#-----------------------------------findPeaks----------------------------------------------
#' Loads the sequence and calculates the coverage of the sequence
#'
#' @param files file,but need index files
#' @param chrs Chromosome information
#' @examples
#' chrs <- c(as.character(1:19),'X','Y')
#' fullcovF <- loadBpCoverages('./dedup_h2.forward.sorted.bam',chrs)
#' @export
loadBpCoverages <- function(files, chrs) {
    chrs <- chrs
    fullCov <- derfinder::fullCoverage(files = files, chrs = chrs)
    file.remove(files)
    return(fullCov)
}


splitPeaks <- function(pos_count, chrom, pos, readlength) {
    chr <- rep(chrom, length(pos_count))
    throshold <- quantile(pos_count, 0.25)
    tempdata <- bumphunter::regionFinder(pos_count, chr = chr, pos = pos, cutoff = throshold, verbose = F)
    tempdata <- subset(tempdata, L > readlength)

    return(tempdata)
}

#' Perform peak calling on a BAM file produced from a scRNA-seq experiment
#'
#' @param fullCov The output of loadBpCoverages
#' @param strand strand
#' @param L Sequence length
#' @param maxwidth Maximum peak width
#' @param cutoff The base-pair level cutoff to use
#' @examples
#' forwardPeaks <-findPeaks(fullcovF,'+',98,1000)
#' reversePeaks <-findPeaks(fullcovR,'-',98,1000)
#' @export
findPeaks <- function(fullCov, Strand, L, maxwidth, cutoff = 0) {
    regionMat <- derfinder::regionMatrix(fullCov = fullCov, L = L, cutoff = cutoff)
    DefPeak <- GenomicRanges::GRanges()
    for (i in names(regionMat)) {
        index <- which(GenomicRanges::width(regionMat[[i]]$regions) > maxwidth)
        tempDataFram <- data.frame()
        for (j in index) {
            rawstart <- start(regionMat[[i]]$regions[j]) - 1
            temp <- regionMat[[i]]$bpCoverage[[j]]
            temp$pos <- c(1:dim(temp)[1])
            first.split.peaks <- splitPeaks(temp$value, i, temp$pos, L)
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
                    cur.split.peaks <- splitPeaks(pos.temp$value, i, pos.temp$pos, L)
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
        DefPeak <- c(DefPeak, regionMat[[i]]$regions[-index])
        if (length(tempDataFram) != 0) {
            DefPeak <- c(DefPeak, regioneR::toGRanges(tempDataFram))
        }
    }
    #DefPeak <- DefPeak[-which(DefPeak$area == GenomicRanges::width(DefPeak))]
    DefPeak <- DefPeak[GenomicRanges::width(DefPeak) > L]
    strand(DefPeak) <- Strand
    if (Strand == "+") {
        DefPeak <- as.data.frame(DefPeak, row.names = 1:length(DefPeak))
        DefPeak <- dplyr::mutate(DefPeak, chr = seqnames, coord = end)
    } else {
        DefPeak <- as.data.frame(DefPeak, row.names = 1:length(DefPeak))
        DefPeak <- dplyr::mutate(DefPeak, chr = seqnames, coord = start)
    }

    return(DefPeak)
}

write.bed <- function(.x, f) {
    write.table(x = .x, file = f, sep = "\t", col.names = F, row.names = F, quote = F)
}

#-----------------------------------GenerateResult/Quantitative----------------------------------------------
#' Generate peak annotation files for quantification
#'
#' @param forwardPeaks The output of \code{\link{findPeaks}},where strand = +
#' @param reversePeaks The output of \code{\link{findPeaks}},where strand = -
#' @param outputdir Output file directory
#' @examples
#' forwardPeaks <-findPeaks(fullcovF,'+',98,1000,cutoff = 0)
#' reversePeaks <-findPeaks(fullcovR,'-',98,1000,cutoff = 0)
#' generateSAF(forwardPeaks,reversePeaks,'./data')
#' @return Generate peaks.saf file in the output directory
#' @export
generateSAF <- function(forwardPeaks, reversePeaks, outputdir) {
    peaks <- rbind(forwardPeaks, reversePeaks)
    peaks$PeakID <- paste0("peak", "_", 1:length(peaks$chr))
    peaks.saf <- peaks[, c("PeakID", "chr", "start", "end", "strand")]

    if (!file.exists(outputdir)) {
        dir.create(outputdir)
    }

    output <- paste0(outputdir, "/", "peaks.saf")
    write.bed(peaks.saf, output)

    return(output)
}

.mergePA <- function(bam, adjBp = 10) {
    peak <- GenomicRanges::GRanges(seq = Rle(bam$chr), ranges = IRanges(bam$coord, bam$coord + adjBp), strand = Rle(bam$strand))

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



#' Generate BAM files for umi_tools to calculate expression
#'
#' @param featureCounts.path The path of the featureCounts
#' @param samtools.path The path of the samtools
#' @param input BAM file
#' @param peakfile The path of the peaks.saf file generated by \code{\link{generateSAF}}
#' @param thread Number of threads
#' @examples
#' input <- './data/demo.bam'
#' peakfile <- './data/peaks.saf'
#' generateFinalBam(featureCounts.path,samtools.path,input,peakfile,24)
#' @return Generate final.bam file in current folder
#' @export
generateFinalBam <- function(featureCounts.path, samtools.path, input, peakfile, thread) {
    output.dir <- dirname(input)
    log <- paste0(output.dir, "/peak_assigned")
    command <- paste0(featureCounts.path, " -a ", peakfile, " -F SAF -t exon -s 1 -M --largestOverlap -o ", log, " -R BAM -T ",
        thread, " ", input)
    system(command = command, wait = T)

    input <- paste0(input, ".featureCounts.bam")
    output <- paste0(dirname(input), "/", "final.bam")
    command <- paste0(samtools.path, " sort -@ ", thread, " ", input, " -o ", output, " && ", samtools.path, " index -@ ",
        thread, " ", output)

    system(command = command, wait = T)

    return(output)
}

#' Calculate the expression level of each peak in each cell
#'
#' @param umitools.path The path of the umi_tools
#' @param input The final.bam file is generated by \code{\link{generateFinalBam}}
#' @param outputdir Output file directory
#' @param TenX Logical value. Is it 10x data ?
#' @examples
#' umitools.path <- '/home/aa/miniconda2/envs/umi_tools/bin/umi_tools'
#' input <- './data/final.bam'
#' outputdir <- './data'
#' countPeaks(umitools.path,input,outputdir,TenX=T)
#' @return Generate counts.tsv.gz file in the output directory
#' @export
countPeaks <- function(umitools.path, input, outputdir, TenX = T) {
    output <- paste0(outputdir, "/counts.tsv.gz")
    if (TenX) {
        command <- paste0(umitools.path, " count --per-gene --gene-tag=XT --assigned-status-tag=XS --method=unique --extract-umi-method=tag --umi-tag=UB --cell-tag=CB --per-cell -I ",
            input, " -S ", output)
        system(command = command, wait = T)
    } else {
        command <- paste0(umitools.path, " count --per-gene --gene-tag=XT --assigned-status-tag=XS --method=unique --per-cell -I ",
            input, " -S ", output)
        system(command = command, wait = T)
    }

    return(output)
}

.findIntersectBypeak <- function(qrypeak, sbjpeak, d = 100) {
    gr1 <- with(qrypeak, GRanges(seqnames = chr, ranges = IRanges(start = coord, end = coord), strand = strand))

    gr2 <- with(sbjpeak, GRanges(seqnames = chr, ranges = IRanges(start = coord, end = coord), strand = strand))

    ov = GenomicRanges::findOverlaps(gr1, gr2, maxgap = d-1, minoverlap = 0L, type = c("any"), select = "all", ignore.strand = FALSE)

    return(qrypeak[unique(ov@from), ])
}

#' Generate single cell expression matrix
#'
#' @param countsfile The decompressed file path of counts.tsv.gz generated by \code{\link{countPeaks}}
#' @param peaksfile The path of the peaks.saf file generated by \code{\link{generateSAF}}
#' @param barcode Vector, and remove the numbers such as -1 in barcode
#' @param tails The result produced by \code{\link{findTails}} or \code{\link{findChrTails}}
#' @param min.cells peak is expressed in at least a few cells,the default value is 2
#' @param min.count peak minimum expression,the default value is 0
#' @examples
#' countsfile <- countPeaks(umitools.path,'./data/final.bam','./data',TenX=T)
#' peakfile <- generateSAF(forwardPeaks,reversePeaks,'./data')
#' tails <- findTails(bamfile = './data/demo.bam')
#' barcode <- read.delim2('./data/barcodes.tsv',header = F)
#' barcode <- gsub('-[0-9]','',barcode$V1)
#' scExpMa <- generatescExpMa(countsfile,peaksfile,barcode,tails,min.cells = 2,min.count = 0)
#' @section Warning:
#' The chromosome information in tails and peaks.saf chromosome information must be consistent.
#' For example, tails chromosome 1 is represented as 1, and peaks.saf chromosome 1 is represented as chr1, this is not possible.
#' it must be 1 or chr1 at the same time.
#' @export
generatescExpMa <- function(countsfile, peaksfile, barcode, tails,d,min.cells = 2,min.count = 0) {
    counts <- read.delim(countsfile, header = T)
    saf <- read.delim(peaksfile, header = F)
    names(saf) <- c("peakID", "chr", "start", "end", "strand")

    counts <- subset(counts, cell %in% barcode)
    peak <- group_by(counts,gene) %>% dplyr::summarize(cells = n(),count = sum(count))
    peak <- subset(peak,cells >= min.cells & count >= min.count)

    saf <- .computePAcoord(saf)
    saf <- .findIntersectBypeak(saf, tails,d=d)

    unionpeak <- intersect(saf$peakID,peak$gene)

    counts <- subset(counts, gene %in% unionpeak)

    saf <- subset(saf,peakID %in% unionpeak)
    scExpMa <- reshape2::dcast(counts, gene ~ cell, value.var = "count")
    scExpMa[is.na(scExpMa)] <- 0


    scExpMa$gene <- factor(scExpMa$gene, levels = saf$peakID)
    scExpMa <- arrange(scExpMa, gene)
    rownames(scExpMa) <- scExpMa$gene
    scExpMa <- scExpMa[, -1]


    scExpMa <- cbind(saf, scExpMa)
    rownames(scExpMa) <- saf$peakID

    return(scExpMa)
}
