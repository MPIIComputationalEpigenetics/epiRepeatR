setClassUnion("BamFileOrNULL", c("BamFile", "NULL"))
setClassUnion("ListOrNULL", c("list", "NULL"))

#' RepeatAlignment Class
#'
#' A general class for storing reads mapped to a reference of repetitive elements
#'
#' @section Slots:
#' \describe{
#'   \item{\code{bamFo}}{
#'       a \code{BamFile} object corresponding to the alignment of reads to a reference
#'       of repetitive elements
#'   }
#'   \item{\code{reference}}{
#'       File name of the reference of REs (fasta file)
#'   }
#'   \item{\code{readCounts}}{
#'       List storing statistics of reads mapped to each reference RE.
#'   } 
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{\code{\link{getReadCounts,RepeatAlignment-method}}}{
#'       Return statistics of reads mapped to each reference RE
#'   }
#'   \item{\code{\link{storeReadCounts,RepeatAlignment-method}}}{
#'       Stored statistics of reads mapped to each reference RE in the object.
#'   }
#'   \item{\code{\link{resetReadCounts,RepeatAlignment-method}}}{
#'       Reset statistics of reads mapped to each reference RE in the object.
#'   }
#' }
#'
#' @name RepeatAlignment-class
#' @rdname RepeatAlignment-class
#' @author Fabian Mueller
#' @noRd
## @exportClass RepeatAlignment
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

#' RepeatAlignment Constructor
#' 
#' @param bamFile	Filename of the BAM file containing the reads aligned to the
#'                  reference of repetitive elements
#' @param reference Filename of the reference of REs (FASTA file)
#' @name RepeatAlignment
#' @rdname RepeatAlignment-class
#' @author Fabian Mueller
#' @noRd
## @export
RepeatAlignment <- function(bamFile, reference=.config$refFasta){
	obj <- new("RepeatAlignment",
		bamFile, reference
	)
	return(obj)
}

if (!isGeneric("getReadCounts")) setGeneric("getReadCounts", function(.obj, ...) standardGeneric("getReadCounts"))
#' getReadCounts-methods
#'
#' Return the sample indentifiers contained in the dataset
#'
#' @param .obj            \code{\linkS4class{RepeatAlignment}} object
#' @param useIdxStats     Optionally specify the the marks which a sample needs to cover in order to be retrieved.
#' @param addGlobalCounts  
#' @return Character vector specifying the sample identifiers
#'
#' @rdname getReadCounts-RepeatAlignment-method
#' @docType methods
#' @aliases getReadCounts
#' @aliases getReadCounts,RepeatAlignment-method
#' @author Fabian Mueller
#' @noRd
## @export
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
					attr(res, "global") <- list(.unmapped=countUnmapped, .mapped=countMapped, .total=countTotal)
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
					attr(res, "global") <- list(.unmapped=countUnmapped, .mapped=countMapped, .total=countTotal)
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
#' storeReadCounts-methods
#'
#' Calls \code{\link{getReadCounts,RepeatAlignment-method}} and stores this information in the \code{\linkS4class{RepeatAlignment}}
#'
#' @param .obj    \code{\linkS4class{RepeatAlignment}} object
#' @param ...     arguments passed to \code{\link{getReadCounts,RepeatAlignment-method}}
#' @return modified \code{\linkS4class{RepeatAlignment}} object
#'
#' @rdname storeReadCounts-RepeatAlignment-method
#' @docType methods
#' @aliases storeReadCounts
#' @aliases storeReadCounts,RepeatAlignment-method
#' @author Fabian Mueller
#' @noRd
## @export
setMethod("storeReadCounts", signature(.obj="RepeatAlignment"),
	function(.obj, ...){
		.obj@readCounts <- getReadCounts(.obj, ...)
		.obj
	}
)
if (!isGeneric("resetReadCounts")) setGeneric("resetReadCounts", function(.obj, ...) standardGeneric("resetReadCounts"))
#' resetReadCounts-methods
#'
#' Reset the read counts in an \code{\linkS4class{RepeatAlignment}} object
#'
#' @param .obj    \code{\linkS4class{RepeatAlignment}} object
#' @return modified \code{\linkS4class{RepeatAlignment}} object
#'
#' @rdname resetReadCounts-RepeatAlignment-method
#' @docType methods
#' @aliases resetReadCounts
#' @aliases resetReadCounts,RepeatAlignment-method
#' @author Fabian Mueller
#' @noRd
## @export
setMethod("resetReadCounts", signature(.obj="RepeatAlignment"),
	function(.obj){
		.obj@readCounts <- NULL
		.obj
	}
)

################################################################################
#' RepeatAlignmentBiSeq Class
#'
#' Daughter class of \code{\linkS4class{RepeatAlignment}} for
#' bisulfite sequencing reads mapped to a reference of repetitive elements
#'
#' @section Slots: see \code{\linkS4class{RepeatAlignment}}
#'
#' @section Methods:
#' \describe{
#'   \item{\code{\link{getMethylationCalls,RepeatAlignmentBiSeq-method}}}{
#'       For each RE in the reference compute the methylation levels in each CpG
#'       as well as read statistics
#'   }
#' }
#'
#' @name RepeatAlignmentBiSeq-class
#' @rdname RepeatAlignmentBiSeq-class
#' @author Fabian Mueller
#' @noRd
## @exportClass RepeatAlignmentBiSeq
setClass("RepeatAlignmentBiSeq",
	contains="RepeatAlignment"
)
#' RepeatAlignmentBiSeq Constructor
#' 
#' @param bamFile	Filename of the BAM file containing the reads aligned to the
#'                  reference of repetitive elements
#' @param reference Filename of the reference of REs (FASTA file)
#' @name RepeatAlignmentBiSeq
#' @rdname RepeatAlignmentBiSeq-class
#' @author Fabian Mueller
#' @noRd
## @export
RepeatAlignmentBiSeq <- function(bamFile, reference=.config$refFasta){
	obj <- new("RepeatAlignmentBiSeq",
		bamFile, reference
	)
	return(obj)
}

if (!isGeneric("getMethylationCalls")) setGeneric("getMethylationCalls", function(.obj, ...) standardGeneric("getMethylationCalls"))
#' getMethylationCalls-methods
#'
#' Get methylation calls for each CpG in each RE along with read statistics
#'
#' @param .obj        \code{\linkS4class{RepeatAlignmentBiSeq}} object
#' @param minBaseQual minumum base quality for a read at the CpG position to be considered in methylation calling
#' @return \code{RepeatMethCalls} object (S3 object). Essentially a list containing a data.frame with methylation calls
#'         (containing position along with methylation, unmethylation and total counts) as well as the total number of reads mapping
#'         to each RE and the number of reads considered in methylation calling
#'
#' @rdname getMethylationCalls-RepeatAlignmentBiSeq-method
#' @docType methods
#' @aliases getMethylationCalls
#' @aliases getMethylationCalls,RepeatAlignmentBiSeq-method
#' @author Fabian Mueller
#' @noRd
## @export
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
################################################################################
#' RepeatAlignmentChip Class
#'
#' Daughter class of \code{\linkS4class{RepeatAlignment}} for
#' reads from enrichment-based experiments mapped to a reference of repetitive elements
#'
#' @section Slots: see \code{\linkS4class{RepeatAlignment}}
#'
#' @section Methods:
#' \describe{
#'   \item{\code{\link{computeEnrichment,RepeatAlignmentChip-method}}}{
#'       Compute the enrichment over Input/WCE for each RE in the reference
#'   }
#' }
#'
#' @name RepeatAlignmentChip-class
#' @rdname RepeatAlignmentChip-class
#' @author Fabian Mueller
#' @noRd
## @exportClass RepeatAlignmentChip
setClass("RepeatAlignmentChip",
	contains="RepeatAlignment"
)
#' RepeatAlignmentChip Constructor
#' 
#' @param bamFile	Filename of the BAM file containing the reads aligned to the
#'                  reference of repetitive elements
#' @param reference Filename of the reference of REs (FASTA file)
#' @name RepeatAlignmentChip
#' @rdname RepeatAlignmentChip-class
#' @author Fabian Mueller
#' @noRd
## @export
RepeatAlignmentChip <- function(bamFile, reference=.config$refFasta){
	obj <- new("RepeatAlignmentChip",
		bamFile,reference
	)
	return(obj)
}
if (!isGeneric("computeEnrichment")) setGeneric("computeEnrichment", function(.obj, inputObj, ...) standardGeneric("computeEnrichment"))
#' computeEnrichment-methods
#'
#' Given another \code{\linkS4class{RepeatAlignmentChip}} object for the Input/WCE, compute the enrichment of the ChIP
#' experiment over the imput for each reference RE
#'
#' @param .obj        \code{\linkS4class{RepeatAlignmentChip}} object
#' @param inputObj    \code{\linkS4class{RepeatAlignmentChip}} object for the Input/WCE
#' @param ...          Arguments passed to \code{storeReadCounts,RepeatAlignment-method}
#' @return \code{RepeatReadEnrichment} object (S3 object). Essentially a list containing the fold change (enrichment) over input,
#'         log2 of that fold change and the number of reads mapping to element in the ChIP and Input experiments
#'         for each RE in the reference
#'
#' @rdname computeEnrichment-RepeatAlignmentChip-method
#' @docType methods
#' @aliases computeEnrichment
#' @aliases computeEnrichment,RepeatAlignmentChip-method
#' @author Fabian Mueller
#' @noRd
## @export
setMethod("computeEnrichment", signature(.obj="RepeatAlignmentChip", inputObj="RepeatAlignment"),
	function(.obj, inputObj, ...){
		eps <- 1e-6
		if (is.null(.obj@readCounts)){
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
			# fc <- (rc.chip[[ss]] * total.input + eps)/(rc.input[[ss]] * total.chip + eps) #overflow problems
			fc <- (rc.chip[[ss]]/total.chip) / (rc.input[[ss]]/total.input)
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
