setClassUnion("ListOrNULL", c("list", "NULL"))
setClassUnion("GenomeRepeatTrackOrNULL", c("GenomeRepeatTrack", "NULL"))

setClass("GenomeAlignment",
	slots = list(
		bamFo = "BamFile",
		repTrack = "GenomeRepeatTrackOrNULL",
		readCounts = "ListOrNULL"
	),
	package = "epiRepeatR"
)

setMethod("initialize", "GenomeAlignment",
	function(.Object,
		bamFile,
		repeatTrack=NULL
	){
		.Object@bamFo <- BamFile(bamFile)
		if (length(index(.Object@bamFo))<1){
			logger.status("Indexing bam file")
			indexBam(.Object@bamFo)
		}
		.Object@repTrack   <- repeatTrack
		.Object@readCounts <- NULL
		.Object
	}
)
GenomeAlignment <- function(bamFile, repeatTrack=NULL){
	object <- new("GenomeAlignment",
		bamFile, repeatTrack
	)
	return(object)
}

if (!isGeneric("getReadCounts")) setGeneric("getReadCounts", function(.obj, ...) standardGeneric("getReadCounts"))
setMethod("getReadCounts", signature(.obj="GenomeAlignment"),
	function(.obj, aggregate=TRUE, addGlobalCounts=TRUE){
		bfo <- .obj@bamFo
		if (is.null(.obj@readCounts)){
			repTrack <- getRepeatInstances(.obj@repTrack)
			if (is.null(repTrack)) logger.error("Need valid repTrack slot when computing read counts for repetitive elements")
			seqNames.bam <- names(scanBamHeader(bfo)$targets)

			seqlvls.stripped <- gsub("^chr", "", seqlevels(repTrack))

			if (sum(seqNames.bam %in% seqlvls.stripped) > sum(seqNames.bam %in% seqlevels(repTrack))){
				logger.info("Stripping 'chr' prefix from reference repeat coordinates")
				seqlevels(repTrack) <- seqlvls.stripped
			}

			excl.seqlvls <- setdiff(seqlevels(repTrack),seqNames.bam)
			if (length(excl.seqlvls)){
				logger.warning(c("The following chroms are not present in the bam file and are disregarded in the analysis:", paste(excl.seqlvls, collapse=",")))
				repTrack <- endoapply(repTrack, FUN=function(x){
					x[!(as.character(seqnames(x)) %in% excl.seqlvls)]
				})
				seqlevels(repTrack) <- setdiff(seqlevels(repTrack), excl.seqlvls)
			}

			repNames <- names(repTrack)
			sbWhat <- scanBamWhat()
			sbFlags <- scanBamFlag(isSecondaryAlignment=FALSE, isNotPassingQualityControls=FALSE)
			readCounts <- lapply(repNames, FUN=function(rn){
				logger.status(c("Counting reads for repeat", rn))
				sbWhich <- repTrack[[rn]]
				# rrs <- scanBam(bfo, param=ScanBamParam(what=sbWhat,which=sbWhich))
				numReads <- countBam(bfo, param=ScanBamParam(flag=sbFlags, what=sbWhat, which=sbWhich))
				res <- sbWhich
				elementMetadata(res) <- data.frame(elementMetadata(res), numReads=numReads$records, nucleotides=numReads$nucleotides)
				return(res)
			})
			names(readCounts) <- repNames
			res <- readCounts
		} else {
			res <- .obj@readCounts
		}
		if (aggregate){
			res <- lapply(res, FUN=function(x){sum(elementMetadata(x)[,"numReads"])})
			if (addGlobalCounts){
				totalC <- countAllReads.bam(Rsamtools::path(bfo))
				mappedEst <- sum(unlist(res))
				# cannot quickly and acurately determine how many reads are mapped or unmapped since one read might overlap multiple instances
				# just provide estimates here
				attr(res, "global") <- list(.unmapped=NA, .mapped=NA, .total=totalC, .mappedEst=mappedEst, .unmappedEst=totalC-mappedEst)
			}
		}
		return(res)
	}
)
if (!isGeneric("storeReadCounts")) setGeneric("storeReadCounts", function(.obj, ...) standardGeneric("storeReadCounts"))
setMethod("storeReadCounts", signature(.obj="GenomeAlignment"),
	function(.obj, ...){
		.obj@readCounts <- getReadCounts(.obj, aggregate=FALSE)
		.obj
	}
)
if (!isGeneric("resetReadCounts")) setGeneric("resetReadCounts", function(.obj, ...) standardGeneric("resetReadCounts"))
setMethod("resetReadCounts", signature(.obj="GenomeAlignment"),
	function(.obj){
		.obj@readCounts <- NULL
		.obj
	}
)

if (!isGeneric("computeEnrichment")) setGeneric("computeEnrichment", function(.obj, inputObj, ...) standardGeneric("computeEnrichment"))
setMethod("computeEnrichment", signature(.obj="GenomeAlignment", inputObj="GenomeAlignment"),
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

if (!isGeneric("normCounts")) setGeneric("normCounts", function(.obj, ...) standardGeneric("normCounts"))
setMethod("normCounts", signature(.obj="GenomeAlignment"),
	function(.obj, method="scale", ...){
		if (is.null(.obj@readCounts)){
			logger.status("Getting read counts")
			.obj <- storeReadCounts(.obj, ...)
		}
		rc <- getReadCounts(.obj)
		rcVec <- unlist(rc)
		totalReads <- sum(rcVec)
		mu <- mean(rcVec, na.rm=TRUE)
		sigma <- sd(rcVec, na.rm=TRUE)

		res <- lapply(names(rc),FUN=function(ss){
			normCount <- NA
			if (method=="scale" && rc[[ss]] > 0){
				normCount <- rc[[ss]]/totalReads
			} else if (method=="zscore" && rc[[ss]] > 0){
				normCount <- (rc[[ss]]-mu)/sigma
			}
			rr <- list(
				normCount=normCount,
				readStats=c(
					numReads=rc[[ss]]
				)
			)	
			return(rr)
		})
		names(res) <- names(rc)
		class(res) <- c("RepeatNormReadCount")
		return(res)
	}
)
