################################################################################
# Loose collection of code that is currently not used in the package
# but might prove useful for validation and in the future
################################################################################

################################################################################
# helper functions to display alignments (so far unused)
################################################################################
expandReadSeqFromCigar <- function(readS,readCigar){
	cigarLen <- as.integer(unlist(strsplit(readCigar,"[MDI]")))
	cigarInstr <- unlist(strsplit(readCigar,"[[:digit:]]+"))
	cigarInstr <- cigarInstr[2:length(cigarInstr)]
	posStart <- 0
	posEnd   <- 0
	res <- ""
	for (i in 1:length(cigarInstr)){
		if (is.element(cigarInstr[i],c("M","I"))){
			posStart <- posEnd + 1
			posEnd   <- posEnd + cigarLen[i]
			res <- paste0(res, substr(readS, posStart, posEnd))
		} else if (is.element(cigarInstr[i],c("D"))){
			res <- paste0(res, paste(rep("-",cigarLen[i]),collapse=""))
		}
	}
	return(res)
}
expandRefSeqFromCigar <- function(refS, mapPos, readCigar){
	cigarLen <- as.integer(unlist(strsplit(readCigar,"[MDI]")))
	cigarInstr <- unlist(strsplit(readCigar,"[[:digit:]]+"))
	cigarInstr <- cigarInstr[2:length(cigarInstr)]
	posStart <- mapPos - 1
	posEnd   <- mapPos - 1
	res <- ""
	for (i in 1:length(cigarInstr)){
		if (is.element(cigarInstr[i],c("M","D"))){
			posStart <- posEnd + 1
			posEnd   <- posEnd + cigarLen[i]
			res <- paste0(res, substr(refS, posStart, posEnd))
		} else if (is.element(cigarInstr[i],c("I"))){
			res <- paste0(res, paste(rep("-",cigarLen[i]),collapse=""))
		}
	}
	return(res)
}
getMatchStr <- function(s1,s2,bisulfite=FALSE){
	if (nchar(s1)!=nchar(s2)) stop("Strings of unequal lengths encountered")
	s1vec <- unlist(strsplit(s1,""))
	s2vec <- unlist(strsplit(s2,""))
	res <- rep(" ",length(s1vec))
	res[s1vec==s2vec] <- "|"
	if (bisulfite){
		res[s1vec=="C" & s2vec=="C"] <- "M"
		res[s1vec=="T" & s2vec=="C"] <- "U"
	}
	res <- paste(res,collapse="")
	return(res)
}
printAln <- function(reads, refSeqs){
	sepLine <- paste(rep("#",80),collapse="")
	subSepLine <- paste(rep("-",80),collapse="")
	for (i in 1:length(reads$qname)){
		readS <- as.character(reads$seq[[i]])
		refS <- as.character(refSeqs[[as.character(reads$rname[i])]])
		readCigar <- reads$cigar[i]
		readSex <- expandReadSeqFromCigar(readS,readCigar)
		refSex  <- expandRefSeqFromCigar(refS, reads$pos[i], readCigar)
		matchS <- getMatchStr(readSex,refSex)
		cat(paste0(
			sepLine,"\n",
			reads$qname[i], " (",as.character(reads$strand[i]),")","\n",
			subSepLine,"\n",
			"READ: ",readSex,"\n",
			"      ",matchS,"\n",
			"REF:  ",refSex,"\n",
			sepLine,"\n"
		))
	}
	invisible(NULL)
}
# #example
# refSeqs <- getDNAStringForReference(.config$refFasta)
# refInfo <- getReferenceInfo(reference)
# refNames <- refInfo[,"id"]
# names(refSeqs) <- refNames
# ss <- "L1MD3"
# alnTargetLengths <- scanBamHeader(bf)$targets
# sbWhat <- scanBamWhat()
# sbWhich <- RangesList(a=IRanges(1,alnTargetLengths[ss]))
# names(sbWhich)[1] <- ss
# rrs <- scanBam(bf, param=ScanBamParam(what=sbWhat,which=sbWhich))
# printAln(rrs[[1]],refSeqs)
