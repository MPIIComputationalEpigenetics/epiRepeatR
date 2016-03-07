#' isCpgAtPosMethInRead
#'
#' check if a single read is methylated for a set of reference CpG positions (coordinates of the C in the reference)
#'
#' @param cPos			vector of cytosine positions
#' @param readS			the read sequence (single character string)
#' @param readQ 		quality scores in phred format. Either a single character string or a vector of integers
#' @param readCigar		the CIGAR mapping string of the read 
#' @param readMapPos	the read's mapping position
#' @param readMapStrand	the read's mapping strand
#' @param minBaseQual  	minimum quality of a base (cytosin and guanine) to be considered for checking CpG methylation status
#' @return a vector of TRUE, FALSE and NA (each element contains the info for a reference position). TRUE if methylated (base is an C), FALSE if unmethylated (base is a T) and NA if quality criteria are not met or the read does not contain a CpG/TpG at that position 
#'
#' @details
#' more details of the individual fields can be found in the \href{https://samtools.github.io/hts-specs/SAMv1.pdf}{SAM/BAM format specifications}
#'
#' @author Fabian Mueller
#' @noRd
#' @examples
#' \donttest{
#' cPos <- c(99,102,106,121,142,146,149,152,200)
#' readMapPos <- 101
#' readMapStrand <- "+"
#' readS <- "XCGXXCAXXXXXXXXXXXXXXXXXXTGXXXXXXXXXXXXTGXXCAXCGXC"
#' readQ <- "DEDDDEDDDDDDDDDDDDDDDDDDDFEDDDDDDDDDDDDEDDDEDDEDDE"
#' readCigar <- "10M5I20M7D15M"
#' # XCGXXCAXXXXXXXXXXXXXXXXXXTGXXXXXXXXXXXXTGXXCAXCGXC
#' # 11111111111111111111111111111111111111111111111111
#' # 00000000011111111111111122222222223334444444444555
#' # 12345678900000012345678901234567890890123456789012
#' isCpgAtPosMethInRead(cPos,readS,readQ,readCigar,readMapPos, readMapStrand, minBaseQual=36)
#' }
isCpgAtPosMethInRead <- function(cPos, readS, readQ, readCigar, readMapPos, readMapStrand, minBaseQual=30){
	if (is.unsorted(cPos)) stop("Sorted vector of positions required")
	if (!is.character(readS) || length(readS) != 1) stop("invalid read sequence argument")
	len <- nchar(readS)
	readSex <- unlist(strsplit(readS,""))

	if (is.character(readQ) && length(readQ) == 1) {
		readQex <- as.integer(PhredQuality(readQ))
	} else if (is.integer(readQ) && length(readQ) == len){
		readQex <- readQ
	} else {
		stop("invalid read quality argument")
	}
	
	#CIGAR parsing
	cigarLen <- as.integer(unlist(strsplit(readCigar,"[MDI]")))
	if (any(is.na(cigarLen))) stop(paste("Inrecognized CIGAR string:",readCigar))
	cigarInstr <- unlist(strsplit(readCigar,"[[:digit:]]+"))
	cigarInstr <- cigarInstr[2:length(cigarInstr)]
	lenCigarCovered <- sum(cigarLen[!is.element(cigarInstr,c("D"))])
	if (lenCigarCovered!=len) stop(paste("CIGAR string (",readCigar,"spans a length of",lenCigarCovered,"while the read sequence is only",len,"long"))

	isMatch <- rep(FALSE,len) #base in the read corresponds to a match to the reference
	numDel <- rep(0,len) # base is not present in read, but only in the reference --> count the number of bases that were deleted before each position in the read
	isIns <- rep(FALSE,len) #base in the read corresponds to an insertion
	posStart <- 0
	posEnd   <- 0
	for (i in 1:length(cigarInstr)){
		if (is.element(cigarInstr[i],c("M"))){
			posStart <- posEnd + 1
			posEnd   <- posEnd + cigarLen[i]
			isMatch[posStart:posEnd] <- TRUE
		} else if (is.element(cigarInstr[i],c("I"))){
			posStart <- posEnd + 1
			posEnd   <- posEnd + cigarLen[i]
			isIns[posStart:posEnd] <- TRUE
		} else if (is.element(cigarInstr[i],c("D"))){
			numDel[posEnd+1] <- cigarLen[i]
		} else {
			stop("Unrecognized CIGAR key")
		}
	}
	cumDel <- cumsum(numDel)
	cumIns <- cumsum(isIns)

	readPosInRef <- 1:len + readMapPos - 1 - cumIns + cumDel

	cPos.read <- match(cPos,readPosInRef)
	gPos.read <- match(cPos+1,readPosInRef)

	cPos.is.cov <- !is.na(cPos.read) & !is.na(gPos.read)

	res <- rep(NA,length(cPos))
	if (readMapStrand=="-"){
		res[cPos.is.cov][readSex[gPos.read[cPos.is.cov]]=="G"] <- TRUE
		res[cPos.is.cov][readSex[gPos.read[cPos.is.cov]]=="A"] <- FALSE
		res[cPos.is.cov][readSex[cPos.read[cPos.is.cov]]!="C"] <- NA
	} else {
		res[cPos.is.cov][readSex[cPos.read[cPos.is.cov]]=="C"] <- TRUE
		res[cPos.is.cov][readSex[cPos.read[cPos.is.cov]]=="T"] <- FALSE
		res[cPos.is.cov][readSex[gPos.read[cPos.is.cov]]!="G"] <- NA 
	}
	#check qual of C and G
	res[cPos.is.cov][readQex[cPos.read[cPos.is.cov]] < minBaseQual] <- NA
	res[cPos.is.cov][readQex[gPos.read[cPos.is.cov]] < minBaseQual] <- NA

	return(res)
	# return(NA)
}
################################################################################
# profiling:
################################################################################
# methCallprofile <- function(){
# 	isMethTab <- do.call("cbind",lapply(1:numReads,FUN=function(i){
# 		isCpgAtPosMethInRead(cPos, as.character(reads$seq[[i]]), as.integer(reads$qual[i]), reads$cigar[i], reads$pos[i], reads$strand[i], minBaseQual=minBaseQual) #careful: only use one pair of [] in reads$qual[i]. otherwise the conversion goes wrong
# 	}))
# }
# methCallprofile.test <- function(N=10000){
# 	i <- 1
# 	minBaseQual <- 30 
# 	cPos <- cPos
	
# 	rseq <- as.character(reads$seq[[i]])
	
# 	rcig <- reads$cigar[i]
# 	rpos <- reads$pos[i]
# 	rstrand <- reads$strand[i]
#	rqual <- as.integer(reads$qual[i])

# 	logger.start("test")
# 	isMethTab <- lapply(1:N,FUN=function(j){
# 		isCpgAtPosMethInRead(cPos, rseq, rqual, rcig, rpos, rstrand, minBaseQual=minBaseQual)
# 	})
# 	logger.completed()
# }
# methCallprofile.test(10000)

# tmp <- tempfile()
# Rprof(tmp, interval = 0.001)
# methCallprofile()
# Rprof(NULL)
# summaryRprof(tmp)
# unlink(tmp)
