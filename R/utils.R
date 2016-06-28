#' getFileTypeFromFileName
#'
#' get the file type from a vector of file names
#'
#' @param fns		vector of filenames
#' @param ignoreZip ignore extions gz,zip,tar,tar.gz und use the main extension
#' @return a vector of file types
#'
#' @author Fabian Mueller
#' @noRd
getFileTypeFromFileName <- function(fns, ignoreZip=FALSE) {
	require(tools)
	if (ignoreZip){
		fns <- gsub("\\.(gz|zip|tar|tar\\.gz)$", "", fns)
	}
	fTypes <- file_ext(fns)
	fTypes[fTypes=="fa"] <- "fasta"
	return(fTypes)
}

#' getCharVecHeadString
#'
#' returns the head of a string vector as comma-separated string
#'
#' @param ss vector of filenames
#' @return a character with the head of the string vector separated by commas
#'
#' @author Fabian Mueller
#' @noRd
getCharVecHeadString <- function(ss, n=5){
	res <- ""
	if (length(ss) < n){
		res <- paste(ss, collapse=", ")
	} else {
		res <- paste(c(ss[1:n], "..."), collapse=", ")
	}
	return(res)
}

#' countAllReads.bam
#'
#' fast count mapped reads to all reference sequence in indexed bamfile using samtools idxstats
#'
#' @param fName		indexed mapped reads file (BAM)
#' @return the number of mapped reads across all reference contigs
#'
#' @author Fabian Mueller
#' @noRd
countAllReads.bam <- function(fName){
	cmd <- paste(.config$samtools.exec,"idxstats",fName)
	res <- NULL
	tryCatch({
			ss <- data.frame(t(data.frame(strsplit(system(cmd, intern=TRUE),"\t"))),stringsAsFactors=FALSE)
			res <- sum(as.numeric(ss[,3])) + sum(as.numeric(ss[,4]))
		},
		error = function(ee) {
			logger.error(c("Not an indexed bam file:",fName,"(",ee$message,")"))
			res <- NULL
		}
	)
	return(res)
}
#' countMappedReads.bam.perSeq
#'
#' fast count the number of mapped reads to each reference contig (and * for unmapped) in indexed bamfile using samtools idxstats
#'
#' @param fName		indexed mapped reads file (BAM)
#' @return a named vector containing the number of reads mapping to each reference contig
#'
#' @author Fabian Mueller
#' @noRd
countMappedReads.bam.perSeq <- function(fName){
	cmd <- paste(.config$samtools.exec,"idxstats",fName)
	res <- NULL
	tryCatch({
			ss <- data.frame(t(data.frame(strsplit(system(cmd, intern=TRUE),"\t"))),stringsAsFactors=FALSE)
			res <- as.numeric(ss[,3])
			names(res) <- ss[,1]
		},
		error = function(ee) {
			logger.error(c("Not an indexed bam file:",fName,"(",ee$message,")"))
			res <- NULL
		}
	)
	return(res)
}
#' countMappedReads.bam
#'
#' fast count mapped reads to all reference sequence (and the * (unmapped) contig) in indexed bamfile using samtools idxstats
#'
#' @param fName		indexed mapped reads file (BAM)
#' @return the sum of all mapped reads across all reference contigs (and the * (unmapped) contig)
#'
#' @author Fabian Mueller
#' @noRd
countMappedReads.bam <- function(fName){
	sum(countMappedReads.bam.perSeq(fName))
}
#' countUnmappedReads.bam
#'
#' fast count unmapped reads in indexed bamfile using samtools idxstats
#'
#' @param fName		indexed mapped reads file (BAM)
#' @return the number of unmapped reads
#'
#' @author Fabian Mueller
#' @noRd
countUnmappedReads.bam <- function(fName){
	cmd <- paste(.config$samtools.exec,"idxstats",fName)
	res <- NULL
	tryCatch({
			ss <- data.frame(t(data.frame(strsplit(system(cmd, intern=TRUE),"\t"))),stringsAsFactors=FALSE)
			res <- sum(as.numeric(ss[,4]))
		},
		error = function(ee) {
			logger.error(c("Not an indexed bam file:",fName,"(",ee$message,")"))
			res <- NULL
		}
	)
	return(res)
}

#' getRepeatMaskerTable
#'
#' obtain the repeat masker track for a given genome from UCSC. NOT USED YET
#'
#' @param species		species. Currently supports only "human" and "mouse"
#' @return UCSC repeat masker table
#'
#' @author Fabian Mueller
#' @noRd
getRepeatMaskerTable <- function(species=.config$species){
	require(rtracklayer)
	session <- browserSession()
	if (species=="human"){
		genome(session) <- "hg19"
		tt <- getTable(ucscTableQuery(session, track="RepeatMasker", table="rmsk"))
	} else if (species=="mouse"){
		genome(session) <- "mm10"
		tt <- getTable(ucscTableQuery(session, track="RepeatMasker", table="rmsk"))
	} else {
		stop("Unknown species")
	}
	return(tt)
}

#' convertPhredCharToInt
#'
#' convert a single character string containing Phred symbols to an integer vector
#'
#' @param ss		single character vector containing Phred symbols (offset:33)
#' @return an integer vector of Phred scores
#'
#' @author Fabian Mueller
#' @noRd
convertPhredCharToInt <- function(ss){
	exploded <- strsplit(ss,"")
	res <- lapply(exploded,FUN=function(x){
		#convert ascii char to integer and subtract 33
		unname(vapply(x,FUN=function(cc){
			strtoi(charToRaw(cc),16L)
		}, integer(1))) - 33L
	})
	return(res)
}

#' colorize.value
#'
#' transform a vector of values to color values
#'
#' @param val		numeric values to be plotted
#' @param rng		output range of colors. Values outside the bounds will be truncated
#' @param colscheme	color panel to be mapped to
#' @return a vector of color values
#'
#' @author Fabian Mueller
#' @noRd
colorize.value <- function(val, rng=c(min(val,na.rm=TRUE),max(val,na.rm=TRUE)), colscheme=colorpanel(100,"white","red")){
    ifelse(val < rng[1], plotval <- rng[1], plotval <- val)
    return(colscheme[round((plotval - rng[1]) / (rng[2] -rng[1])  * (length(colscheme)-1),0)+1])
}

#' getFgColorForBg
#'
#' get the appropriate foreground color for given background colors. "white" for dark background, "black" for light background
#'
#' @param bg		 a vector of background color values
#' @param thresh	mean rgb color intensity threshold for switching colors
#' @return a vector of foreground color values
#'
#' @author Fabian Mueller
#' @noRd
getFgColorForBg <- function(bg,thresh=100){
    col.rgb.vec <- as.vector(col2rgb(bg))
    ifelse(mean(col.rgb.vec)<thresh,return("white"),return("black"))
}

#' bitwCompl
#'
#' get the bitwise complement of an integer representation of a number, i.e. negate
#' the bitwise representation
#'
#' @param a an integer representation of a binary
#' @return integer representation of the bitwise complement
#'
#' @author Fabian Mueller
#' @noRd
bitwCompl <- function(a){
	2^ceiling(log2(a))-a-1
}

#' explainSamFlags
#'
#' Convert SAM flag integer(s) into bit representation
#'
#' @param x  (vector of) integer representation(s) of SAM flag(s)
#' @return a bit matrix containing explanations
#'
#' @author Fabian Mueller
#' @noRd
explainSamFlags <- function(x){
	xbit <- do.call("rbind", lapply(x, FUN=function(xc) {as.logical(intToBits(xc))}))
	xbit <- xbit[,1:12,drop=FALSE]
	colnames(xbit) <- c(
		"read_paired",
		"read_properPair",
		"read_unmapped",
		"mate_unmapped",
		"read_reverse",
		"mate_reverse",
		"read_read1_inpair",
		"read_read2_inpair",
		"read_secondary_alignment",
		"read_fails_qc",
		"read_duplicate",
		"read_supplementary_alignment"
	)
	return(xbit)
}

#' getFlagStats
#'
#' gets the samtools flagstat results for a given bam file.
#'
#' @param bamFn bam file
#' @return ...
#'
#' @author Fabian Mueller
#' @noRd
getFlagStats <- function(bamFn){
	require(stringr)
	selFlags <- c(
		"total"      = "in total",
		"mapped"     = "mapped",
		"paired"     = "paired in sequencing",
		"properPair" = "properly paired",
		"read1"      = "read1",
		"read2"      = "read2",
		"duplicate"  = "duplicates",
		"secondary"  = "secondary",
		"singletons" = "singletons"
	)
	parseFlagStats <- function(ss){
		ss <- str_trim(gsub("\\(.+\\)", "", ss))
		sl <- strsplit(ss," ")
		tt <- do.call("rbind", lapply(sl, FUN=function(x){
			if (x[2] != "+") stop(paste0("error in pearsing flagstats for file ",bamFn, ": second element in line is non '+'"))
			return(c(x[1], x[3], paste(x[4:length(x)], collapse=" ")))
		}))
		tt <- data.frame(
			qc_pass=as.integer(tt[,1]),
			qc_fail=as.integer(tt[,2]),
			description=tt[,3],
			stringsAsFactors=FALSE
		)
		tt$total <- tt$qc_pass + tt$qc_fail
		return(tt)
	}
	cmd <- paste(.config$samtools.exec,"flagstat",bamFn)
	res <- NULL
	tryCatch({
			res <- parseFlagStats(system(cmd, intern=TRUE))
		},
		error = function(ee) {
			logger.error(c("Error parsing flagstats for file:", bamFn, "(",ee$message,")"))
			res <- NULL
		}
	)
	selflagInDesc <- selFlags %in% res$description
	if (!all(selflagInDesc)){
		stop(paste0("The following flags were not found in the flagstats for file", bamFn, paste(names(selFlags)[!selflagInDesc],collapse=",")))
	}
	res <- res[!duplicated(res$description),c("qc_pass", "qc_fail", "total", "description")]
	rownames(res) <- res$description
	res <- res[selFlags,]
	rownames(res) <- names(selFlags)
	return(res)
}

#' getReadStatsFromSample
#'
#' sample a bam file and get read statistics from that sample
#'
#' @param bamFn bam file
#' @param frac  fraction of reads to be sampled
#' @return a list containing read statistics
#'
#' @author Fabian Mueller
#' @noRd
getReadStatsFromSample <- function(bamFn, frac=0.001){
	# tmpBamFn <- tempfile(pattern="sampledBam_", fileext=".bam")
	ss <- system2(.config$samtools.exec, c("view", "-s", frac, bamFn), stdout=TRUE)
	ss <- strsplit(ss, "\t")
	readLens <- vapply(ss, FUN=function(x){nchar(x[10])}, integer(1))
	flags <- vapply(ss, FUN=function(x){as.integer(x[2])}, integer(1))
	flagTab <- explainSamFlags(flags)
	flagRates <- colSums(flagTab)/nrow(flagTab) 
	res <- list(
		readLength.summary = summary(readLens),
		flagRates = flagRates
	)
	return(res)
}


#' plotRepeatAlignmentStats
#'
#' plot the repeat alignment statistic from the table output by the pipeline step "repeatAlignmentStats"
#'
#' @param x  filename for the table or data frame containing alignment statistics
#' @return nothing of particular intest
#'
#' @author Fabian Mueller
#' @export
plotRepeatAlignmentStats <- function(x){
	if (is.character(x)){
		if (file.exists(x)){
			tt <- read.table(x, sep="\t", header=TRUE, stringsAsFactors=FALSE)
		} else {
			stop("invalid alignment stats file")
		}
	} else if (is.data.frame(x)){
		tt <- x
	} else {
		stop("invalid table argument")
	}
	df2p <- tt[,c("sampleName", "mark", "dataType", "totalReads", "mappedReads", "mappingRate")]
	pp <- ggplot(df2p) + aes(x=totalReads, y=mappingRate, color=mark, shape=dataType) + geom_point() #+ geom_text(aes(label=sampleName))	
	pp
}

#' mergeMethGr
#'
#' Merge methylation and total counts from both strands and identical coordinates
#'
#' @param gr GRanges object with methylation and total read counts annotated with T and M in the metadata
#' @return a GRanges object containing summed methylation calls from both strands for identical coordinates
#'
#' @author Fabian Mueller
#' @noRd
mergeMethGr <- function(gr){
	gr <- sort(gr, ignore.strand=TRUE)
	# djb <- disjointBins(gr, ignore.strand=TRUE)
	# create GRanges object with unique coordinates
	gr.dj <- gr
	strand(gr.dj) <- "*"
	elementMetadata(gr.dj) <- NULL
	gr.dj <- unique(gr.dj)
	oo <- findOverlaps(gr.dj, gr, type="equal", ignore.strand=TRUE)
	tt <- elementMetadata(gr)
	tt <- tt[subjectHits(oo),]
	sumsM <- as.integer(tapply(tt$M, queryHits(oo), FUN=sum))
	sumsT <- as.integer(tapply(tt$T, queryHits(oo), FUN=sum))
	elementMetadata(gr.dj) <- data.frame(M=sumsM, T=sumsT)
	return(gr.dj)
}
#' parseMcTable.bissnp
#'
#' Parse methylation calls from BisSNP output format
#'
#' @param mcFile methylation call file
#' @return a GRanges object containing methylation calls for cytosine positions
#'
#' @author Fabian Mueller
#' @noRd
parseMcTable.bissnp <- function(mcFile){
	tt <- read.table(mcFile, sep="\t", header=FALSE, skip=1)
	colnames(tt)[1:6] <- c("chrom", "start", "end", "perc", "T", "strand")
	doShift <- tt$strand != "-"
	tt$start[doShift] <- tt$start[doShift] + 1L
	tt$end[doShift] <- tt$end[doShift] + 1L

	gr <- GRanges(tt$chrom, IRanges(tt$start, width=2), strand=tt$strand)
	elementMetadata(gr) <- data.frame(M=as.integer(round(tt$perc/100 * tt$T)),T=tt$T)

	gr.merged <- mergeMethGr(gr)
	return(gr.merged)
}

#' logger.cmd.args
#'
#' Log command line options
#'
#' @param cmdArgs command line options as returned by \code{parse_args} of the \code{argparse} package
#' @return nothing of particular interest
#'
#' @author Fabian Mueller
#' @export
logger.cmd.args <- function(cmdArgs){
	require(stringr)
	argNames <- names(cmdArgs)
	argStrings <- as.character(cmdArgs)
	argStrings <- paste0(str_pad(argNames, max(nchar(argNames)), side="right"), " : ", argStrings)
	logger.start("Parameter settings")
		for (i in 1:length(argStrings)){
			logger.info(argStrings[i])
		}
	logger.completed()
	invisible(NULL)
}

#' normalizeMatrix
#'
#' Apply normalizatio to columns of a given matrix
#'
#' @param x      a matrix of values to be normalized
#' @param method method of normalization. Currently supported are:
#'               \code{none} (no normalization),
#'               \code{standard} (subtract the mean, devide by standard deviation),
#'               \code{scale} (scale to the interval [0,1]) and
#'               \code{quantile} (Quantile normalization)
#' @param ...    arguments passed down to the actual normalization functions
#' @return the normalized matrix
#'
#' @author Fabian Mueller
#' @export
normalizeMatrix <- function(x, method="standard", ...){
	require(matrixStats)
	normFuns <- list(
		none = function(x){
			return(x)
		},
		standard = function(x){
			rr <- t( (t(x) - colMeans(x, na.rm=TRUE))/colSds(x, na.rm=TRUE) )
			return(rr)
		},
		quantile = function(x){
			require(preprocessCore)
			# xr <- apply(x, 2, FUN=function(cc){rank(cc, ties.method="min")})
			# xs <- apply(x, 2, sort)
			# xo <- apply(x, 2, FUN=function(cc){order(cc)})
			rr <- normalize.quantiles(x)
			return(rr)
		},
		scale = function(x, a=0, b=1){
			rr <- a + t( ((t(x) - colMins(x, na.rm=TRUE))*(b-a))/(colMaxs(x, na.rm=TRUE)-colMins(x, na.rm=TRUE)) )
			return(rr)
		}
	)
	if (!is.matrix(x)) logger.error("expected matrix [normalizeMatrix]")
	if (!is.element(method, names(normFuns))) logger.error(c("unknown normalization method:", method))
	res <- normFuns[[method]](x, ...)
	return(res)
}