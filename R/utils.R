#' getFileTypeFromFileName
#'
#' get the file type from a vector of file names
#'
#' @param fns		vector of filenames
#' @return a vector of file types
#'
#' @author Fabian Mueller
#' @noRd
getFileTypeFromFileName <- function(fns) {
	require(tools)
	fTypes <- file_ext(fns)
	fTypes[fTypes=="fa"] <- "fasta"
	return(fTypes)
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
