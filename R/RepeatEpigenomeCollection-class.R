setClassUnion("ListOrNULL", c("list", "NULL"))

setClass("RepeatEpigenomeCollection",
	slots = list(
		annot="data.frame",
		repRef="RepeatReference",
		samples="character",
		marks="character",
		markTypes="factor",
		epiQuant="list"
	),
	package = "epiRepeatR"
)


#' inferMarkTypes
#'
#' given a character vector of epigenetic marks, returns a factor containing
#' data types belonging to that mark or raises an error if mark is unknown
#'
#' @param marks character vector of mark names
#' @return factor of data types
#'
#' @author Fabian Mueller
#' @noRd
inferMarkTypes <- function(marks){
	patternHistoneChip <- "^H[34][K][0-9]+(me|ac)[0-9]?$"
	markTypes <- ifelse(
		marks=="DNAmeth",
		"DNAmeth",
		ifelse(
			grepl(patternHistoneChip, marks, ignore.case=TRUE),
			"ChIPseq",
			NA
		)
	)
	if (any(is.na(markTypes))){
		stop(paste0("Unable to infer mark type for mark(s):", paste(marks[is.na(markTypes)], collapse=",")))
	}
	return(factor(markTypes))
}

setMethod("initialize", "RepeatEpigenomeCollection",
	function(.Object,
		quantFns,
		sampleNames,
		markNames,
		annot,
		repRef=RepeatReference()
	){
		if (length(quantFns)!=length(sampleNames) || length(quantFns)!=length(markNames)) {
			stop("incompatible quantFns, sampleNames and markNames parameters. Should all be equal in length")
		}
		if (is.null(rownames(annot))){
			stop("Invalid annotation table: no rownames")
		}
		samples <- rownames(annot)
		if (!all(sampleNames %in% samples)){
			stop(paste0(
				"The following sample names are not present in the annotation table: ",
				paste(setdiff(unique(sampleNames),samples), collapse=",")
			))
		}
		marks <- sort(unique(markNames))
		markTypes <- inferMarkTypes(marks)

		epiQuant <- lapply(samples, FUN=function(sn){
			resList <- lapply(marks, FUN=function(mn){
				res <- NULL
				fn <- quantFns[sampleNames==sn & markNames==mn]
				if (length(fn) > 0){
					if(length(fn) > 1){
						stop(paste0("Multiple quantification file names specified for ", sn, " -- ", mn))
					}
					res <- readRDS(fn[1])
				}
				return(res)
			})
			names(resList) <- marks
			return(resList)
		})
		names(epiQuant) <- samples

		.Object@samples   <- samples
		.Object@marks     <- marks
		.Object@markTypes <- markTypes
		.Object@annot     <- annot
		.Object@repRef    <- repRef
		.Object@epiQuant  <- epiQuant
		
		.Object
	}
)

RepeatEpigenomeCollection <- function(quantFns, sampleNames, markNames, annot, repRef=RepeatReference()){
	obj <- new("RepeatEpigenomeCollection",
		quantFns, sampleNames, markNames, annot, repRef
	)
	return(obj)
}

setMethod("show", "RepeatEpigenomeCollection",
	function(object) {
		sns <- getSamples(object)
		mns <- getMarks(object)
		rr <- getRepRef(object)

		cat("RepeatEpigenomeCollection\n")
		cat(c("#samples : ", length(sns), " [", getCharVecHeadString(sns), "]", "\n"), sep="")
		cat(c("#marks   : ", length(mns), " [", getCharVecHeadString(mns), "]", "\n"), sep="")
		cat(c("#repeats : ", length(getRepeatIds(rr)), "\n"), sep="")
	}
)

if (!isGeneric("getSamples")) setGeneric("getSamples", function(.Object) standardGeneric("getSamples"))
setMethod("getSamples", signature(.Object="RepeatEpigenomeCollection"),
	function(.Object){
		return(.Object@samples)
	}
)
if (!isGeneric("getMarks")) setGeneric("getMarks", function(.Object) standardGeneric("getMarks"))
setMethod("getMarks", signature(.Object="RepeatEpigenomeCollection"),
	function(.Object){
		return(.Object@marks)
	}
)
if (!isGeneric("getAnnot")) setGeneric("getAnnot", function(.Object) standardGeneric("getAnnot"))
setMethod("getAnnot", signature(.Object="RepeatEpigenomeCollection"),
	function(.Object){
		return(.Object@annot)
	}
)
if (!isGeneric("getRepRef")) setGeneric("getRepRef", function(.Object) standardGeneric("getRepRef"))
setMethod("getRepRef", signature(.Object="RepeatEpigenomeCollection"),
	function(.Object){
		return(.Object@repRef)
	}
)

if (!isGeneric("getRepeatScores")) setGeneric("getRepeatScores", function(.Object, ...) standardGeneric("getRepeatScores"))
setMethod("getRepeatScores", signature(.Object="RepeatEpigenomeCollection"),
	function(.Object, mark, dropEmptySamples=FALSE){
		if (!is.element(mark, getMarks(.Object))){
			stop(paste0("unknown mark:", mark))
		}
		markType <- inferMarkTypes(mark)
		scoreFun <- function(x){NA}
		if (markType == "DNAmeth"){
			scoreFun <- function(x){
				mean(x$methCalls[,"numM"]/x$methCalls[,"numT"], na.rm=TRUE)
			}
		} else if (markType == "ChIPseq"){
			scoreFun <- function(x){
				x$log2fc
			}
		} else {
			stop(paste0("Unknown data type for mark:",mark))
		}

		repRefNames <- getRepeatIds(getRepRef(.Object))
		nReps <- length(repRefNames)
		scoreMat <- do.call("cbind",lapply(getSamples(.Object), FUN=function(sn){
			sres <- rep(NA, nReps)
			quantRes <- .Object@epiQuant[[sn]][[mark]]
			if (!is.null(quantRes)){
				sres <- sapply(quantRes[repRefNames], FUN=function(r){
					if (is.null(r)) return(NA)
					return(scoreFun(r))
				})
			}
			return(sres)
		}))
		colnames(scoreMat) <- getSamples(.Object)
		rownames(scoreMat) <- repRefNames
		if (dropEmptySamples){
			scoreMat <- scoreMat[, !colAlls(is.na(scoreMat))]
		}
		return(scoreMat)
	}
)

if (!isGeneric("getRepeatCovg")) setGeneric("getRepeatCovg", function(.Object, ...) standardGeneric("getRepeatCovg"))
setMethod("getRepeatCovg", signature(.Object="RepeatEpigenomeCollection"),
	function(.Object, mark, dropEmptySamples=FALSE){
		if (!is.element(mark, getMarks(.Object))){
			stop(paste0("unknown mark:", mark))
		}
		markType <- inferMarkTypes(mark)
		covgFun <- function(x){NA}
		if (markType == "DNAmeth"){
			covgFun <- function(x){
				x$readStats["numReads"]
			}
		} else if (markType == "ChIPseq"){
			covgFun <- function(x){
				x$readStats["numReads_chip"]
			}
		} else {
			stop(paste0("Unknown data type for mark:",mark))
		}

		repRefNames <- getRepeatIds(getRepRef(.Object))
		nReps <- length(repRefNames)
		covgMat <- do.call("cbind",lapply(getSamples(.Object), FUN=function(sn){
			sres <- rep(NA, nReps)
			quantRes <- .Object@epiQuant[[sn]][[mark]]
			if (!is.null(quantRes)){
				sres <- sapply(quantRes[repRefNames], FUN=function(r){
					if (is.null(r)) return(NA)
					return(covgFun(r))
				})
			}
			return(sres)
		}))
		colnames(covgMat) <- getSamples(.Object)
		rownames(covgMat) <- repRefNames
		if (dropEmptySamples){
			covgMat <- covgMat[, !colAlls(is.na(covgMat))]
		}
		return(covgMat)
	}
)
#NOTE: getRepeatCovg and getRepeatScores can lead to different repeat-sample combinations to be NA
# for methylation, not all reads covering the repeat might have CpG information
# for ChIPseq the coverage might be 0 for either the input or the chip

################################################################################
# Sandbox
################################################################################
if (FALSE){
anaDir <- "/DEEP_fhgfs/projects/fmueller/repeatEpigenetics/epiRepeatR/analysis/deepBlood_v02mergedInput"
am <- readRDS(file.path(anaDir, "config", "anaMan.rds"))

inFileTable <- read.table(file.path(anaDir, "config", "plotRepeatMarkTree_inputFiles.tsv"), sep="\t", comment.char="", header=TRUE, stringsAsFactors=FALSE)

quantFns <- inFileTable[,"fileName"]
sampleNames <- inFileTable[,"sampleName"]
markNames <- inFileTable[,"markName"]
annot <- getSampleAnnot(am)

rec <- RepeatEpigenomeCollection(quantFns, sampleNames, markNames, annot)
# mark <- "H3K9me3"
mark <- "DNAmeth"
sm <- getRepeatScores(rec, mark)
cm <- getRepeatCovg(rec, mark)
}
