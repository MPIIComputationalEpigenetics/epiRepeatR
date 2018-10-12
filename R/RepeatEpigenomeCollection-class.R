setClassUnion("ListOrNULL", c("list", "NULL"))

#' RepeatEpigenomeCollection Class
#'
#' A class for storing quantitative epigenetic data in repetitive elements of the genome
#' 
#' @details
#' Multiple epigenetic marks are stored and each repetitive element (RE) in a given reference is assigned quantitative
#' scores for each mark. A dataset of multiple samples can be annotated.
#'
#' @section Slots:
#' \describe{
#'   \item{\code{annot}}{
#'       Annotation for the samples in the dataset.
#'   }
#'   \item{\code{repRef}}{
#'       Repeat reference (of class \code{RepeatReference})
#'   }
#'   \item{\code{samples}}{
#'       Samples contained in the dataset.
#'   }
#'   \item{\code{marks}}{
#'       epigenetic marks contained in the dataset.
#'   }
#'   \item{\code{markTypes}}{
#'       Annotation for the samples in the dataset.
#'   }
#'   \item{\code{epiQuant}}{
#'       Data structure for storing the quantifications of the epigenetic marks across marks and samples
#'   }
#' 
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{\code{\link{getSamples,RepeatEpigenomeCollection-method}}}{
#'       Return the sample indentifiers contained in the dataset.
#'   }
#'   \item{\code{\link{getMarks,RepeatEpigenomeCollection-method}}}{
#'       Return the epigenetic marks contained in the dataset.
#'   }
#'   \item{\code{\link{getAnnot,RepeatEpigenomeCollection-method}}}{
#'       Return the sample sample annotation table.
#'   }
#'   \item{\code{\link{getSampleMarkTable,RepeatEpigenomeCollection-method}}}{
#'       Return a table containing sample-mark combinations covered in the dataset.
#'   }
#'   \item{\code{\link{getRepRef,RepeatEpigenomeCollection-method}}}{
#'       Retrieve values of epigenetic quantifications for each RE in each sample.
#'   }
#'   \item{\code{\link{getRepeatScores,RepeatEpigenomeCollection-method}}}{
#'       Return the sample indentifiers contained in the dataset.
#'   }
#'   \item{\code{\link{getRepeatCovg,RepeatEpigenomeCollection-method}}}{
#'       Retrieve the number of reads covering each RE in each sample.
#'   }
#'   \item{\code{\link{filterRepRefMeth,RepeatEpigenomeCollection-method}}}{
#'       Filter RE that do not fulfill certain coverage criteria in all of the contained
#'       bisulfite datasets.
#'   }
#'   \item{\code{\link{filterRepRefChip,RepeatEpigenomeCollection-method}}}{
#'       Filter RE that do not fulfill certain coverage criteria in all of the contained
#'       enrichment datasets.
#'   }
#' }
#'
#' @name RepeatEpigenomeCollection-class
#' @rdname RepeatEpigenomeCollection-class
#' @author Fabian Mueller
#' @exportClass RepeatEpigenomeCollection
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
	markTypes <- marks
	markTypes[grepl(patternHistoneChip, markTypes)] <- "ChIPseq"
	markTypes[markTypes=="accessibility"] <- "Acc"
	markTypes[markTypes=="ATACseq"] <- "Acc"
	isValid <- markTypes %in% c("DNAmeth", "ChIPseq", "Acc")
	if (any(!isValid)){
		stop(paste0("Unable to infer mark type for mark(s):", paste(unique(marks[!isValid]), collapse=",")))
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

#' RepeatEpigenomeCollection Constructor
#' 
#' @param quantFns	      Vector of filenames which contains quantification results
#' @name RepeatEpigenomeCollection
#' @rdname RepeatEpigenomeCollection-class
#' @author Fabian Mueller
#' @export
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

################################################################################
# Direct Accessors
################################################################################

if (!isGeneric("getSamples")) setGeneric("getSamples", function(.Object, ...) standardGeneric("getSamples"))
#' getSamples-methods
#'
#' Return the sample indentifiers contained in the dataset
#'
#' @param .Object \code{\linkS4class{RepeatEpigenomeCollection}} object
#' @param marks   Optionally specify the the marks which a sample needs to cover in order to be retrieved.
#' @return Character vector specifying the sample identifiers
#'
#' @rdname getSamples-RepeatEpigenomeCollection-method
#' @docType methods
#' @aliases getSamples
#' @aliases getSamples,RepeatEpigenomeCollection-method
#' @author Fabian Mueller
#' @export
setMethod("getSamples", signature(.Object="RepeatEpigenomeCollection"),
	function(.Object, marks=NULL){
		sns <- .Object@samples
		if (is.null(marks)){
			return(sns)
		} else {
			smt <- getSampleMarkTable(.Object)
			if (!all(marks %in% getMarks(.Object))){
				logger.warning(c("RepeatEpigenomeCollection object does not contain marks:", paste(setdiff(marks, getMarks(.Object)), collapse=","), "--> no samples returned"))
				return(c())
			}
			snInds <- rowAlls(smt[,marks, drop=FALSE])
			return(sns[snInds])
		}
	}
)
if (!isGeneric("getMarks")) setGeneric("getMarks", function(.Object) standardGeneric("getMarks"))
#' getMarks-methods
#'
#' Return the epigenetic marks contained in the dataset
#'
#' @param .Object \code{\linkS4class{RepeatEpigenomeCollection}} object
#' @return Character vector specifying the epigenetic marks contained in the dataset
#'
#' @rdname getMarks-RepeatEpigenomeCollection-method
#' @docType methods
#' @aliases getMarks
#' @aliases getMarks,RepeatEpigenomeCollection-method
#' @author Fabian Mueller
#' @export
setMethod("getMarks", signature(.Object="RepeatEpigenomeCollection"),
	function(.Object){
		return(.Object@marks)
	}
)
if (!isGeneric("getAnnot")) setGeneric("getAnnot", function(.Object) standardGeneric("getAnnot"))
#' getAnnot-methods
#'
#' Return the sample sample annotation table
#'
#' @param .Object \code{\linkS4class{RepeatEpigenomeCollection}} object
#' @return Sample annotation as \code{data.frame}
#'
#' @rdname getAnnot-RepeatEpigenomeCollection-method
#' @docType methods
#' @aliases getAnnot
#' @aliases getAnnot,RepeatEpigenomeCollection-method
#' @author Fabian Mueller
#' @export
setMethod("getAnnot", signature(.Object="RepeatEpigenomeCollection"),
	function(.Object){
		return(.Object@annot)
	}
)
if (!isGeneric("getSampleMarkTable")) setGeneric("getSampleMarkTable", function(.Object) standardGeneric("getSampleMarkTable"))
#' getSampleMarkTable-methods
#'
#' Return a table containing sample-mark combinations covered in the dataset
#'
#' @param .Object \code{\linkS4class{RepeatEpigenomeCollection}} object
#' @return Logical matrix of dimension numSamples X numMarks indicating whether a given sample-mark combination is covered
#'
#' @rdname getSampleMarkTable-RepeatEpigenomeCollection-method
#' @docType methods
#' @aliases getSampleMarkTable
#' @aliases getSampleMarkTable,RepeatEpigenomeCollection-method
#' @author Fabian Mueller
#' @export
setMethod("getSampleMarkTable", signature(.Object="RepeatEpigenomeCollection"),
	function(.Object){
		sns <- .Object@samples
		marks <- getMarks(.Object)
		hasMark <- do.call("rbind", lapply(sns, FUN=function(sn){
			sapply(marks, FUN=function(mn){
				!is.null(.Object@epiQuant[[sn]][[mn]])
			})
		}))
		rownames(hasMark) <- sns
		return(hasMark)
	}
)
if (!isGeneric("getRepRef")) setGeneric("getRepRef", function(.Object) standardGeneric("getRepRef"))
#' getRepRef-methods
#'
#' Return the reference repeat object
#'
#' @param .Object \code{\linkS4class{RepeatEpigenomeCollection}} object
#' @return The repeat reference as \code{\linkS4class{RepeatReference}} object
#'
#' @rdname getRepRef-RepeatEpigenomeCollection-method
#' @docType methods
#' @aliases getRepRef
#' @aliases getRepRef,RepeatEpigenomeCollection-method
#' @author Fabian Mueller
#' @export
setMethod("getRepRef", signature(.Object="RepeatEpigenomeCollection"),
	function(.Object){
		return(.Object@repRef)
	}
)

################################################################################
# Getting Repeat Scores and Coverage
################################################################################

if (!isGeneric("getRepeatScores")) setGeneric("getRepeatScores", function(.Object, ...) standardGeneric("getRepeatScores"))
#' getRepeatScores-methods
#'
#' Retrieve values of epigenetic quantifications for each RE in each sample
#'
#' @param .Object          \code{\linkS4class{RepeatEpigenomeCollection}} object
#' @param mark             Epigenetic mark for which the score information is retrieved
#' @param dropEmptySamples Logical indicating whether samples in which the mark is not present should be dropped
#' @param minCpGcov        [\code{DNAmeth} only] number of reads covering each CpG in order for it to be considered in the
#'                         mean methylation level
#' @return A matrix of dimension numRE X numSamples containing scores. Scores are mean methylation levels for DNA methylation and 
#'         log2 of the fold change for enrichment based methods.
#'
#' @rdname getRepeatScores-RepeatEpigenomeCollection-method
#' @docType methods
#' @aliases getRepeatScores
#' @aliases getRepeatScores,RepeatEpigenomeCollection-method
#' @author Fabian Mueller
#' @export
setMethod("getRepeatScores", signature(.Object="RepeatEpigenomeCollection"),
	function(.Object, mark, dropEmptySamples=FALSE, minCpGcov=getConfigElement("meth.minCpGcov")){
		if (!is.element(mark, getMarks(.Object))){
			stop(paste0("unknown mark:", mark))
		}
		markType <- inferMarkTypes(mark)
		scoreFun <- function(x){NA}
		if (markType == "DNAmeth"){
			if (is.null(minCpGcov) || minCpGcov < 1){
				scoreFun <- function(x){
					mean(x$methCalls[,"numM"]/x$methCalls[,"numT"], na.rm=TRUE)
				}
			} else {
				scoreFun <- function(x){
					hasCov <- x$methCalls[,"numT"] >= minCpGcov
					if (sum(hasCov) > 0){
						return(mean(x$methCalls[hasCov,"numM"]/x$methCalls[hasCov,"numT"], na.rm=TRUE))
					} else {
						return(NA)
					}
				}
			}
		} else if (markType == "ChIPseq"){
			scoreFun <- function(x){
				x$log2fc
			}
		} else if (markType == "Acc"){
			scoreFun <- function(x){
				x$normCount
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
#' getRepeatCovg-methods
#'
#' Retrieve the number of reads covering each RE in each sample
#'
#' @param .Object          \code{\linkS4class{RepeatEpigenomeCollection}} object
#' @param mark             Epigenetic mark for which the coverage information is retrieved
#' @param dropEmptySamples Logical indicating whether samples in which the mark is not present should be dropped
#' @param type             Coverage type. Should be one of the following:
#'                         \code{"numReads"}: [default] number of reads covering a RE
#'                         \code{"maxCpGcov"}: [\code{DNAmeth} only] maximum number of reads covering a single CpG
#'                         \code{"numInstances"}: [\code{DNAmeth} from genome methylation calls only] number of repeat instances in the genome
#' @return A matrix of dimension numRE X numSamples containing coverage values
#'
#' @rdname getRepeatCovg-RepeatEpigenomeCollection-method
#' @docType methods
#' @aliases getRepeatCovg
#' @aliases getRepeatCovg,RepeatEpigenomeCollection-method
#' @author Fabian Mueller
#' @export
setMethod("getRepeatCovg", signature(.Object="RepeatEpigenomeCollection"),
	function(.Object, mark, dropEmptySamples=FALSE, type="numReads"){
		if (!is.element(mark, getMarks(.Object))){
			stop(paste0("unknown mark:", mark))
		}
		markType <- inferMarkTypes(mark)
		covgFun <- function(x){NA}
		if (markType == "DNAmeth"){
			if (type=="maxCpGcov"){
				covgFun <- function(x){
					if (nrow(x$methCalls) < 1) return(NA)
					return(max(x$methCalls$numT, na.rm=TRUE))
				}
			} else if (type=="numInstances"){
				# for data inferred from genome methylation calls only: use the
				# number of unique repeat instances
				covgFun <- function(x){
					if (is.element("repeatInstanceIndex", colnames(x$methCalls))){
						return(length(unique(x$methCalls$repeatInstanceIndex)))
					} else {
						return(0)
					}
					
				}
			} else if (type=="numReads") {
				covgFun <- function(x){
					x$readStats["numReads"]
				}
			} else {
				logger.error("Invalid coverage type")
			}
		} else if (markType == "ChIPseq"){
			covgFun <- function(x){
				x$readStats["numReads_chip"]
			}
		} else if (markType == "Acc"){
			covgFun <- function(x){
				x$readStats["numReads"]
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
# Differential analysis
################################################################################
if (!isGeneric("getComparisonInfo")) setGeneric("getComparisonInfo", function(.Object, ...) standardGeneric("getComparisonInfo"))
#' getComparisonInfo-methods
#'
#' Retrieve values of differential epigenetic quantifications for each RE in each sample
#'
#' @param .Object          \code{\linkS4class{RepeatEpigenomeCollection}} object
#' @param annotCol         name of the column in the sample annotation table on which the comparison is based
#' @param name.grp1        name of the first group in the corresponding column of the sample annotation table
#' @param name.grp2        name of the second group in the corresponding column of the sample annotation table
#' @param cmpName          (optional) strings specifying a name for the comparison
#' @param covariateCols    column names in the sample annotation table of potentially confounding covariates for the analysis
#' @return Comparison information as structured \code{S3} object of class \code{comparisonInfo}
#'
#' @rdname getComparisonInfo-RepeatEpigenomeCollection-method
#' @docType methods
#' @aliases getComparisonInfo
#' @aliases getComparisonInfo,RepeatEpigenomeCollection-method
#' @author Fabian Mueller
#' @export
setMethod("getComparisonInfo", signature(.Object="RepeatEpigenomeCollection"),
	function(.Object, annotCol, name.grp1, name.grp2, cmpName=NULL, covariateCols=character(0)){
		annotTab <- getAnnot(.Object)
		if (!is.element(annotCol, colnames(annotTab))) stop(paste0("Unknown sample annotation column:", annotCol))
		if (!all(covariateCols %in% colnames(annotTab))) stop(paste0("Not all covariate columns are contained in the sample annotation table"))

		covariateCols <- setdiff(covariateCols, annotCol)
		if (is.null(cmpName)){
			cmpName <- paste(annotCol, "cmp", name.grp1, "vs", name.grp2, sep="_")
		}
		sampleIdx.grp1 <- which(annotTab[,annotCol]==name.grp1)
		if (length(sampleIdx.grp1) < 1) stop(paste0("No samples found for group ", name.grp1, " in sample annotation column ", annotCol))
		sampleIdx.grp2 <- which(annotTab[,annotCol]==name.grp2)
		if (length(sampleIdx.grp2) < 1) stop(paste0("No samples found for group ", name.grp2, " in sample annotation column ", annotCol))

		cmpInfo <- list(
			cmpName=cmpName,
			name.grp1=name.grp1,
			name.grp2=name.grp2,
			sampleIdx.grp1=sampleIdx.grp1,
			sampleIdx.grp2=sampleIdx.grp2,
			annotCol=annotCol,
			covariateCols=covariateCols,
			designF=as.formula(paste0("~", paste(c(covariateCols, annotCol),collapse="+")))
		)
		class(cmpInfo) <- "comparisonInfo"
		return(cmpInfo)
	}
)

if (!isGeneric("getComparisons")) setGeneric("getComparisons", function(.Object, ...) standardGeneric("getComparisons"))
#' getComparisons-methods
#'
#' Retrieve a list of \code{comparisonInfo} objects based
#'
#' @param .Object          \code{\linkS4class{RepeatEpigenomeCollection}} object
#' @param annotCols        names of the column in the sample annotation table on which the comparisons are based.
#'                         \code{NULL} (default) indicates that all columns will be inspected for compatible comparisons
#' @param minGroupSize     minimum number of samples to be considered a group
#' @param maxGroupCount    maximum number of groups in each sample annotation column
#' @param allPairwise      if more than two groups exist for a sample annotation column, should all pairwise comparisons be conducted (as opposed to one vs all).
#'                         currently the only allowed value is \code{TRUE} (default)
#' @param covariateCols    column names in the sample annotation table of potentially confounding covariates for the analysis
#' @return List of comparison information objects (structured \code{S3} object of class \code{comparisonInfo})
#'
#' @rdname getComparisons-RepeatEpigenomeCollection-method
#' @docType methods
#' @aliases getComparisons
#' @aliases getComparisons,RepeatEpigenomeCollection-method
#' @author Fabian Mueller
#' @export
setMethod("getComparisons", signature(.Object="RepeatEpigenomeCollection"),
	function(.Object, annotCols=NULL, minGroupSize=2L, maxGroupCount=2, allPairwise=TRUE, covariateCols=character(0)){
		annotTab <- getAnnot(.Object)
		if (length(annotCols) < 1) annotCols <- colnames(annotTab)
		res <- list()
		for (cc in annotCols){
			grpCounts <- table(annotTab[,cc])
			grpCounts <- grpCounts[grpCounts>minGroupSize]
			nGroups <- length(grpCounts)
			if (nGroups > 1 && nGroups <= maxGroupCount){
				if (nGroups == 2){
					res <- c(res, list(getComparisonInfo(.Object, cc, names(grpCounts)[1], names(grpCounts)[2], covariateCols=covariateCols)))
				} else {
					if (allPairwise){
						groupNamePairs <- t(combn(names(grpCounts),2))
						for (i in 1:nrow(groupNamePairs)){
							res <- c(res, list(getComparisonInfo(.Object, cc, groupNamePairs[i,1], groupNamePairs[i,2], covariateCols=covariateCols)))
						}
					} else {
						stop("one vs all comparisons are not yet supported")
					}
				}
			}
		}
		if (length(res) > 0) {
			names(res) <- sapply(res, FUN=function(x){x$cmpName})
		}
		return(res)
	}
)

if (!isGeneric("getRepeatScoresDiff")) setGeneric("getRepeatScoresDiff", function(.Object, ...) standardGeneric("getRepeatScoresDiff"))
#' getRepeatScoresDiff-methods
#'
#' Retrieve values of differential epigenetic quantifications for each RE in each sample
#'
#' @param .Object          \code{\linkS4class{RepeatEpigenomeCollection}} object
#' @param mark             Epigenetic mark for which the score information is retrieved
#' @param compInfo         \code{S3} object of class \code{comparisonInfo} as returned by \code{\link{getComparisonInfo,RepeatEpigenomeCollection-method}}
#'                         or \code{\link{getComparisons,RepeatEpigenomeCollection-method}}
#' @return A matrix with differential statistics for each repeat subfamily.
#'
#' @rdname getRepeatScoresDiff-RepeatEpigenomeCollection-method
#' @docType methods
#' @aliases getRepeatScoresDiff
#' @aliases getRepeatScoresDiff,RepeatEpigenomeCollection-method
#' @author Fabian Mueller
#' @export
setMethod("getRepeatScoresDiff", signature(.Object="RepeatEpigenomeCollection"),
	function(.Object, mark, compInfo){
		if (!is.element(mark, getMarks(.Object))){
			stop(paste0("unknown mark:", mark))
		}
		if (class(compInfo) != "comparisonInfo") {
			stop("Invalid compInfo paramater. Expected 'comparisonInfo' object")
		}
		markType <- inferMarkTypes(mark)
		scoreFun <- function(x){NA}
		if (markType == "DNAmeth"){
			stop(paste0("Don't know how to compute differential methylation yet"))
		} else if (is.element(markType, c("ChIPseq", "Acc"))){
			countMat <- getRepeatCovg(.Object, mark)
			dm <- computeDiffMat.counts.deseq2(countMat, getAnnot(.Object), compInfo$designF, compInfo$annotCol, compInfo$name.grp1, compInfo$name.grp2)
			diffScores <- data.frame(
				diffScore=dm[,"log2FoldChange"],
				diffPval=dm[,"pvalue"],
				diffPval.adj=dm[,"padj"],
				diffRank=dm[,"cRank_rerank"],
				score.g1=log2(dm[,"baseMean.g1"]),
				score.g2=log2(dm[,"baseMean.g2"])
			)
		} else {
			stop(paste0("Unknown data type for mark:",mark))
		}
		return(diffScores)
	}
)

################################################################################
# Filtering
################################################################################

if (!isGeneric("filterRepRefMeth")) setGeneric("filterRepRefMeth", function(.Object, ...) standardGeneric("filterRepRefMeth"))
#' filterRepRefMeth-methods
#'
#' Given a \code{\linkS4class{RepeatEpigenomeCollection}} object, remove repeats that do not fulfill
#' the coverage criteria in any sample
#'
#' @param .Object	      \code{\linkS4class{RepeatEpigenomeCollection}} object
#' @param minReads 		  threshold for the minimum number of reads required to cover a repeat
#' @param minCpGs		  threshold for the minimum number of CpGs that must be contained in a repeat
#' @param minCpGcov		  threshold specifying that mininum number of reads that must cover at least one
#'                        CpG in the repeat in order for it to be retained 
#' @param covgType        type of read coverage to be applied (important for filtering via \code{minReads}).
#'                        See \code{\link{getRepeatCovg,RepeatEpigenomeCollection-method}} for possible values.
#' @return modified \code{\linkS4class{RepeatEpigenomeCollection}} object
#' 
#' @details
#' A repeat must fulfill the criteria in all samples in order to be retained
#'
#' @rdname filterRepRefMeth-RepeatEpigenomeCollection-method
#' @docType methods
#' @aliases filterRepRefMeth
#' @aliases filterRepRefMeth,RepeatEpigenomeCollection-method
#' @author Fabian Mueller
#' @export
setMethod("filterRepRefMeth", signature(.Object="RepeatEpigenomeCollection"),
	function(.Object,  minReads=getConfigElement("plotRepTree.meth.minReads"), minCpGs=getConfigElement("plotRepTree.meth.minCpGs"), minCpGcov=getConfigElement("meth.minCpGcov"), covgType="numReads"){
		res <- .Object
		repRef <- getRepRef(.Object)
		repRefNames <- getRepeatIds(repRef)
		nReps <- length(repRefNames)
		zeroVec <- rep(0, nReps)
		methSamples <- getSamples(.Object, marks="DNAmeth")
		if (length(methSamples) < 1){
			logger.warning("No sample with DNA methylation measurements found. --> skipping filtering")
			return(.Object)
		}
		survive.numCG <- matrix(TRUE,nrow=length(repRefNames),ncol=length(methSamples))
		if (minCpGs > 1){
			numCG <- do.call("cbind",lapply(.Object@epiQuant[methSamples], FUN=function(x){
				rr <- zeroVec
				if (length(x[["DNAmeth"]][repRefNames]) > 0) {
					rr <- sapply(x[["DNAmeth"]][repRefNames], FUN=function(r){
						if (is.null(r)) return(0) else return(length(r$methCalls$cPos))
					})
				}
				return(rr)
			}))
			survive.numCG <- numCG >= minCpGs
		}

		survive.numReads <- matrix(TRUE,nrow=length(repRefNames),ncol=length(methSamples))
		if (minReads > 1){
			# numReads <- do.call("cbind",lapply(.Object@epiQuant[methSamples], FUN=function(x){
			# 	rr <- zeroVec
			# 	if (length(x[["DNAmeth"]][repRefNames]) > 0) {
			# 		rr <- sapply(x[["DNAmeth"]][repRefNames], FUN=function(r){
			# 			if (is.null(r)) return(0) else return(r$readStats["numReads"])
			# 		})
			# 	}
			# 	return(rr)
			# }))
			numReads <- getRepeatCovg(.Object, "DNAmeth", type=covgType)
			if (sum(numReads>0, na.rm=TRUE)){
				survive.numReads <- numReads >= minReads
			} else {
				logger.warning("No element has a number of associated reads > 0. Ignore this warning if you are working with methylation calls originating from genome alignments.")
			}
		}

		survive.CpGcov <- matrix(TRUE,nrow=length(repRefNames),ncol=length(methSamples))
		if (!(is.null(minCpGcov) || minCpGcov < 1)){
			maxCovg <- do.call("cbind",lapply(.Object@epiQuant[methSamples], FUN=function(x){
				rr <- zeroVec
				if (length(x[["DNAmeth"]][repRefNames]) > 0) {
					rr <- sapply(x[["DNAmeth"]][repRefNames], FUN=function(r){
						if (is.null(r) || all(is.na(r$methCalls$numT))) return(0) else return(max(r$methCalls$numT, na.rm=TRUE))
					})
				}
				return(rr)
			}))
			survive.CpGcov <- maxCovg >= minCpGcov
		}
		survive <- survive.numCG & survive.numReads & survive.CpGcov
		survive <- apply(survive,1,all)
		repRef.filtered <- filterRepeats_wl(repRef, repRefNames[survive])

		res@repRef <- repRef.filtered
		return(res)
	}
)

if (!isGeneric("filterRepRefChip")) setGeneric("filterRepRefChip", function(.Object, ...) standardGeneric("filterRepRefChip"))
#' filterRepRefChip-methods
#'
#' Given a \code{\linkS4class{RepeatEpigenomeCollection}} object, remove repeats that do not fulfill
#' the coverage criteria in any sample
#'
#' @param .Object	      \code{\linkS4class{RepeatEpigenomeCollection}} object
#' @param minReads 		  threshold for the minimum number of reads required to cover a repeat
#' @return modified \code{\linkS4class{RepeatEpigenomeCollection}} object
#' 
#' @details
#' A repeat must fulfill the criteria in all samples in order to be retained
#'
#' @rdname filterRepRefChip-RepeatEpigenomeCollection-method
#' @docType methods
#' @aliases filterRepRefChip
#' @aliases filterRepRefChip,RepeatEpigenomeCollection-method
#' @author Fabian Mueller
#' @export
setMethod("filterRepRefChip", signature(.Object="RepeatEpigenomeCollection"),
	function(.Object,  minReads=getConfigElement("plotRepTree.meth.minReads")){
		res <- .Object
		repRef <- getRepRef(.Object)
		repRefNames <- getRepeatIds(repRef)
		nReps <- length(repRefNames)
		zeroVec <- rep(0, nReps)

		chipMarks <- getMarks(.Object)[inferMarkTypes(getMarks(.Object))=="ChIPseq"]
		if (length(chipMarks) < 1){
			logger.warning("No ChIP mark found")
			return(res)
		}
		smt <- getSampleMarkTable(.Object)[, chipMarks, drop=FALSE]
		chipSamples <- getSamples(.Object)[rowAnys(smt)]
		if (length(chipSamples) < 1){
			logger.warning("No sample with ChIP measurements found. --> skipping filtering")
			return(.Object)
		}
		numReadsTabs.chip <- lapply(chipMarks, FUN=function(mn){
			do.call("cbind",lapply(.Object@epiQuant[chipSamples], FUN=function(x){
				rr <- zeroVec
				if (length(x[[mn]][repRefNames]) > 0) {
					rr <- sapply(x[[mn]][repRefNames], FUN=function(r){
						if (is.null(r)) return(0) else return(r$readStats["numReads_chip"])
					})
				}
				return(rr)
			}))
		})
		names(numReadsTabs.chip) <- chipMarks
		numReadsTabs.input <- lapply(chipMarks, FUN=function(mn){
			do.call("cbind",lapply(.Object@epiQuant[chipSamples], FUN=function(x){
				rr <- zeroVec
				if (length(x[[mn]][repRefNames]) > 0) {
					rr <- sapply(x[[mn]][repRefNames], FUN=function(r){
						if (is.null(r)) return(0) else return(r$readStats["numReads_input"])
					})
				}
				return(rr)
			}))
		})
		names(numReadsTabs.input) <- chipMarks

		# keep repeats in which for any chip mark, all samples fullfill the read coverage criterion for both chip and input
		survive <- matrix(FALSE, nrow=nReps, ncol=length(chipSamples))
		for (mn in chipMarks){
			survive <- survive | (numReadsTabs.chip[[mn]] >= minReads & numReadsTabs.input[[mn]] >= minReads)
		}
		survive <- apply(survive,1,all)
		repRef.filtered <- filterRepeats_wl(repRef, repRefNames[survive])

		res@repRef <- repRef.filtered
		return(res)
	}
)

if (!isGeneric("filterRepRefAcc")) setGeneric("filterRepRefAcc", function(.Object, ...) standardGeneric("filterRepRefAcc"))
#' filterRepRefAcc-methods
#'
#' Given a \code{\linkS4class{RepeatEpigenomeCollection}} object, remove repeats that do not fulfill
#' the coverage criteria in any sample
#'
#' @param .Object	      \code{\linkS4class{RepeatEpigenomeCollection}} object
#' @param minReads 		  threshold for the minimum number of reads required to cover a repeat
#' @return modified \code{\linkS4class{RepeatEpigenomeCollection}} object
#' 
#' @details
#' A repeat must fulfill the criteria in all samples in order to be retained
#'
#' @rdname filterRepRefAcc-RepeatEpigenomeCollection-method
#' @docType methods
#' @aliases filterRepRefAcc
#' @aliases filterRepRefAcc,RepeatEpigenomeCollection-method
#' @author Fabian Mueller
#' @export
setMethod("filterRepRefAcc", signature(.Object="RepeatEpigenomeCollection"),
	function(.Object,  minReads=getConfigElement("plotRepTree.meth.minReads")){
		res <- .Object
		repRef <- getRepRef(.Object)
		repRefNames <- getRepeatIds(repRef)
		nReps <- length(repRefNames)
		zeroVec <- rep(0, nReps)

		accMarks <- getMarks(.Object)[inferMarkTypes(getMarks(.Object))=="Acc"]
		if (length(accMarks) < 1){
			logger.warning("No accessibility mark found")
			return(res)
		}
		smt <- getSampleMarkTable(.Object)[, accMarks, drop=FALSE]
		accSamples <- getSamples(.Object)[rowAnys(smt)]
		if (length(accSamples) < 1){
			logger.warning("No sample with accessibility measurements found. --> skipping filtering")
			return(.Object)
		}
		numReadsTabs.acc <- lapply(accMarks, FUN=function(mn){
			do.call("cbind",lapply(.Object@epiQuant[accSamples], FUN=function(x){
				rr <- zeroVec
				if (length(x[[mn]][repRefNames]) > 0) {
					rr <- sapply(x[[mn]][repRefNames], FUN=function(r){
						if (is.null(r)) return(0) else return(r$readStats["numReads"])
					})
				}
				return(rr)
			}))
		})
		names(numReadsTabs.acc) <- accMarks

		# keep repeats in which for any acc mark, all samples fullfill the read coverage criterion
		survive <- matrix(FALSE, nrow=nReps, ncol=length(accSamples))
		for (mn in accMarks){
			survive <- survive | numReadsTabs.acc[[mn]] >= minReads
		}
		survive <- apply(survive,1,all)
		# keep repeats that have NA scores for at least one sample in any acc mark
		survive.na <- matrix(FALSE, nrow=nReps, ncol=length(accSamples))
		for (mn in accMarks){
			scs <- getRepeatScores(.Object, mn)
			survive.na <- survive.na | !rowAlls(is.na(scs))
		}
		survive <- survive & survive.na

		repRef.filtered <- filterRepeats_wl(repRef, repRefNames[survive])

		res@repRef <- repRef.filtered
		return(res)
	}
)
