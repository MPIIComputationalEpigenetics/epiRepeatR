setClassUnion("BamFileOrNULL", c("BamFile", "NULL"))
setClassUnion("ListOrNULL", c("list", "NULL"))

setClass("RepeatAlignment",
	slots = list(
		bamFo = "BamFile",
		reference = "character",
		readCounts = "ListOrNULL"
	),
	package = "epiRepeatR"
)

setMethod("initialize", "RepeatAlignment",
	function(.Object,
		bamFile,
		reference=.config$refFasta
	){
		.Object@bamFo <- BamFile(bamFile)
		if (length(index(.Object@bamFo))<1){
			logger.status("Indexing bam file")
			indexBam(.Object@bamFo)
		}
		.Object@reference <- reference
		.Object@readCounts <- NULL
		.Object
	}
)
RepeatAlignment <- function(bamFile, reference=.config$refFasta){
	obj <- new("RepeatAlignment",
		bamFile, reference
	)
	return(obj)
}

if (!isGeneric("getReadCounts")) setGeneric("getReadCounts", function(.obj, ...) standardGeneric("getReadCounts"))
setMethod("getReadCounts", signature(.obj="RepeatAlignment"),
	function(.obj, useIdxStats=FALSE, addGlobalCounts=TRUE){
		if (is.null(.obj@readCounts)){
			bf <- .obj@bamFo
			refInfo <- getReferenceInfo(.obj@reference)
			refNames <- refInfo[,"id"]
			alnTargetLengths <- scanBamHeader(bf)$targets
			alnNames <- names(alnTargetLengths)
			if (length(setdiff(alnNames,refNames)) > 0){
				logger.error("There are repeats in the alignment that are not present in the reference")
			}
			if (useIdxStats){
				#faster version using bam idxstats
				bamFn <- Rsamtools::path(bf)
				mappedReadCounts <- countMappedReads.bam.perSeq(bamFn)
				res <- lapply(refNames,FUN=function(ss){
					if (is.element(ss,names(mappedReadCounts))){
						curCount <- unname(mappedReadCounts[ss])
					} else {
						curCount <- 0L
					}
					return(curCount)
				})
				names(res) <- refNames
				if (addGlobalCounts){
					countMapped <- sum(mappedReadCounts)
					countUnmapped <- countUnmappedReads.bam(bamFn)
					countTotal <- countMapped + countUnmapped
					res <- c(res, list(.unmapped=countUnmapped, .mapped=countMapped, .total=countTotal))
				}
			} else {
				#more elaborate version allowing for better filtering
				#prepare scanbam parameters
				sbWhat <- scanBamWhat()
				res <- lapply(refNames,FUN=function(ss){
					logger.status(c("Processing sequence:",ss))
					# print(ss)
					sbWhich <- RangesList(a=IRanges(1,alnTargetLengths[ss]))
					names(sbWhich)[1] <- ss
					curCount <- countBam(bf, param=ScanBamParam(what=sbWhat,which=sbWhich))$records
					return(curCount)
				})
				names(res) <- refNames
				if (addGlobalCounts){
					# countTotal    <- countBam(bf, param=ScanBamParam(what=sbWhat))$records
					countMapped   <- countBam(bf, param=ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE), what=sbWhat))$records
					countUnmapped <- countBam(bf, param=ScanBamParam(flag=scanBamFlag(isUnmappedQuery=TRUE),  what=sbWhat))$records
					countTotal <- countMapped + countUnmapped
					res <- c(res, list(.unmapped=countUnmapped, .mapped=countMapped, .total=countTotal))
				}
			}
			# class(res) <- c("RepeatReadCounts")
		} else {
			res <- .obj@readCounts
		}
		return(res)
	}
)
if (!isGeneric("storeReadCounts")) setGeneric("storeReadCounts", function(.obj, ...) standardGeneric("storeReadCounts"))
setMethod("storeReadCounts", signature(.obj="RepeatAlignment"),
	function(.obj, ...){
		.obj@readCounts <- getReadCounts(.obj, ...)
		.obj
	}
)
if (!isGeneric("resetReadCounts")) setGeneric("resetReadCounts", function(.obj, ...) standardGeneric("resetReadCounts"))
setMethod("resetReadCounts", signature(.obj="RepeatAlignment"),
	function(.obj){
		.obj@readCounts <- NULL
		.obj
	}
)

setClass("RepeatAlignmentBiSeq",
	contains="RepeatAlignment"
)
RepeatAlignmentBiSeq <- function(bamFile, reference=.config$refFasta){
	obj <- new("RepeatAlignmentBiSeq",
		bamFile, reference
	)
	return(obj)
}

if (!isGeneric("getMethylationCalls")) setGeneric("getMethylationCalls", function(.obj, ...) standardGeneric("getMethylationCalls"))
setMethod("getMethylationCalls", signature(.obj="RepeatAlignmentBiSeq"),
	function(.obj, minBaseQual=30){
		bf <- .obj@bamFo
		refSeqs <- getDNAStringForReference(.obj@reference)
		refInfo <- getReferenceInfo(.obj@reference)
		refNames <- refInfo[,"id"]
		names(refSeqs) <- refNames
		refCgSites <- vmatchPattern("CG",refSeqs)
		refCgStarts <- lapply(refCgSites,FUN=start)

		alnTargetLengths <- scanBamHeader(bf)$targets
		alnNames <- names(alnTargetLengths)
		if (length(setdiff(alnNames,refNames)) > 0){
			logger.error("There are repeats in the alignment that are not present in the reference")
		}
		#prepare scanbam parameters
		sbWhat <- scanBamWhat()
		#prepare pileup parameters
		# puParam <- PileupParam(max_depth=.Machine$integer.max, include_deletions=FALSE, min_base_quality=minBaseQual)#, distinguish_strands=FALSE)
		res <- lapply(refNames,FUN=function(ss){
			logger.status(c("Processing sequence:",ss))
			# print(ss)
			cPos <- refCgStarts[[ss]]
			if (length(cPos)<1) return(NULL)
			sbWhich <- RangesList(a=IRanges(1,alnTargetLengths[ss])) #list of coordinates for reads to be extract
			names(sbWhich)[1] <- ss
			rrs <- scanBam(bf, param=ScanBamParam(what=sbWhat,which=sbWhich))
			reads <- rrs[[1]] #only one chromosome/repeat included in rrs
			numReads <- length(reads$qname)
			# posRes <- lapply(1:length(cPos),FUN=function(i){
			# 	pos <- cPos[i]
			# 	sbWhich <- RangesList(a=IRanges(pos,pos+1)) #list of coordinates for reads to be extract
			# 	names(sbWhich)[1] <- ss
			# 	sbParam <- ScanBamParam(what=sbWhat,which=sbWhich)
			# 	pu <- pileup(bf, scanBamParam=sbParam, pileupParam=puParam)
			# 	resi <- pu
			# 	return(resi)
			# })
			if (numReads<1) return(NULL)
			rseqs <- as.character(reads$seq)
			# rquals.org <- sapply(1:numReads,FUN=function(i){as.integer(reads$qual[i])}) #takes up much time
			rquals <- convertPhredCharToInt(as.character(reads$qual)) #speed up
			rcigars <- reads$cigar
			rposs <- reads$pos
			rstrands <- reads$strand
			# # REMARK: due to overhead, parallel is even slower than non-parallel: don't parallelize over reads here
			# if (parallel.isEnabled() && numReads>1){
			# 	mbq <- minBaseQual #for passing this parameter on to the parallel environment, it (strangely) has to be assigned explicitely
			# 	# logger.start("par")
			# 	isMethTab <- foreach(i=1:numReads, .combine='cbind',.multicombine=TRUE,.maxcombine=200) %dopar% {
			# 		epiRepeatR:::isCpgAtPosMethInRead(cPos, rseqs[i], rquals[[i]], rcigars[i], rposs[i], rstrands[i], minBaseQual=mbq)
			# 	}
			# 	# logger.completed() #.export=c("epiRepeatR:::isCpgAtPosMethInRead","reads","minBaseQual","numReads")
			# } else {
				# logger.start("nopar")
				isMethTab <- do.call("cbind",lapply(1:numReads,FUN=function(i){
					isCpgAtPosMethInRead(cPos, rseqs[i], rquals[[i]], rcigars[i], rposs[i], rstrands[i], minBaseQual=minBaseQual)
				}))
				# logger.completed()
			# }
			numM <- rowSums(isMethTab, na.rm=TRUE)
			numU <- rowSums(!isMethTab, na.rm=TRUE)
			numT <- numM + numU
			isReadInformative <- colSums(is.na(isMethTab)) < length(cPos)
			methCalls <- data.frame(
				cPos=cPos,
				numM=numM,
				numU=numU,
				numT=numT
			)
			logger.info(c("found", sum(isReadInformative), "informative reads (of", numReads, "reads aligned). Called methylation for", sum(numT > 0), "of", length(numT), "CpGs [", ss, "]"))
			ress <- list(
				methCalls=methCalls,
				readStats=c("numReads"=numReads, "numReads_informative"=sum(isReadInformative))
			)
			return(ress)
		})
		names(res) <- refNames
		class(res) <- c("RepeatMethCalls")
		return(res)
	}
)
#example:
#ss <- "LTR2752"

setClass("RepeatAlignmentChip",
	contains="RepeatAlignment"
)
RepeatAlignmentChip <- function(bamFile, reference=.config$refFasta){
	obj <- new("RepeatAlignmentChip",
		bamFile,reference
	)
	return(obj)
}
if (!isGeneric("computeEnrichment")) setGeneric("computeEnrichment", function(.obj, inputObj, ...) standardGeneric("computeEnrichment"))
setMethod("computeEnrichment", signature(.obj="RepeatAlignmentChip", inputObj="RepeatAlignment"),
	function(.obj, inputObj, ...){
		eps <- 1e-6
		if (is.null(inputObj@readCounts)){
			logger.status("getting read counts for ChIP")
			.obj <- storeReadCounts(.obj, ...)
		}
		if (is.null(inputObj@readCounts)){
			logger.status("getting read counts for input")
			inputObj <- storeReadCounts(inputObj, ...)
		}
		rc.chip <- getReadCounts(.obj)
		rc.input <- getReadCounts(inputObj)
		total.chip <- sum(unlist(rc.chip))
		total.input <- sum(unlist(rc.input))
		if (!setequal(names(rc.chip), names(rc.input))) stop("Different sequences represented in input and ChIP objects")
		res <- lapply(names(rc.chip),FUN=function(ss){
			fc <- (rc.chip[[ss]] * total.input + eps)/(rc.input[[ss]] * total.chip + eps)
			if (rc.input[[ss]]==0) fc <- NA
			rr <- list(
				fc=fc,
				log2fc=log2(fc),
				readStats=c(
					numReads_chip=rc.chip[[ss]],
					numReads_input=rc.input[[ss]]
				)
			)	
			return(rr)
		})
		names(res) <- names(rc.chip)
		class(res) <- c("RepeatReadEnrichment")
		return(res)
	}
)

#helper functions to display alignment
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

getAlnStats <- function(bamFn.aln, bamFn.unaln){
	ra <- RepeatAlignment(bamFn.aln)
	rcl.aln <- getReadCounts(ra, useIdxStats=TRUE, addGlobalCounts=TRUE)

	rc.mapped <- rcl.aln[[".mapped"]]
	rc.total <- countAllReads.bam(bamFn.unaln)
	res <- data.frame(
		totalReads=rc.total,
		mappedReads=rc.mapped,
		mappingRate=rc.mapped/rc.total
	)
	return(res)
}

# #example
# refSeqs <- getDNAStringForReference(.config$refFasta)
# refInfo <- getReferenceInfo(reference)
# refNames <- refInfo[,"id"]
# names(refSeqs) <- refNames
# ss <- "L1MD3"
# sbWhat <- scanBamWhat()
# sbWhich <- RangesList(a=IRanges(1,alnTargetLengths[ss]))
# names(sbWhich)[1] <- ss
# rrs <- scanBam(bf, param=ScanBamParam(what=sbWhat,which=sbWhich))
# printAln(rrs[[1]],refSeqs)
