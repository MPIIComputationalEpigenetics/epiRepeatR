setClassUnion("ListOrNULL", c("list", "NULL"))

#' RepeatReference Class
#'
#' A class for storing annotation on the reference of repetitive elements (REs)
#' 
#' @section Slots:
#' \describe{
#'   \item{\code{reference}}{
#'       File name of the reference of REs (fasta file)
#'   }
#'   \item{\code{repeatInfo}}{
#'       data.frame containing annotation (columns) for each RE (rows).
#'       This information is typically derived from the annotation in the
#'       reference FASTA file.
#'   }
#'   \item{\code{repeatInfoList}}{
#'       List od additional information for each RE. This information is
#'       typically obtained from the annotaiton contained in RepBaseUpdate's
#'       EMBL files.
#'   }
#'   \item{\code{sequences}}{
#'       \code{DNAStringSet} containing sequence information for each RE.
#'   }
#'   \item{\code{sequences}}{
#'       Species the reference pertains to (as string)
#'   }
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{\code{\link{getRepeatInfo,RepeatReference-method}}}{
#'       Blubb
#'   }
#' }
#'
#' @name RepeatReference-class
#' @rdname RepeatReference-class
#' @author Fabian Mueller
#' @noRd
## @exportClass RepeatReference
setClass("RepeatReference",
	slots = list(
		reference="character",
		repeatInfo="data.frame",
		repeatInfoList="ListOrNULL",
		sequences="DNAStringSet",
		species="character"
	),
	package = "epiRepeatR"
)

setMethod("initialize", "RepeatReference",
	function(.Object,
		reference=.config$refFasta,
		species=.config$species
	){
		.Object@reference=reference
		.Object@repeatInfo=getReferenceInfo(reference)
		.Object@repeatInfoList=NULL
		.Object@sequences=getDNAStringForReference(reference)
		.Object@species=species
		.Object
	}
)

#' RepeatReference Constructor
#' 
#' @param reference Filename of the reference of REs (FASTA file)
#' @param species   Species the reference pertains to
#' @name RepeatReference
#' @rdname RepeatReference-class
#' @author Fabian Mueller
#' @noRd
## @export
RepeatReference <- function(reference=.config$refFasta, species=.config$species){
	obj <- new("RepeatReference",
		reference
	)
	obj <- addRepeatInfoFromEmbl(obj)
	return(obj)
}

if (!isGeneric("addRepeatInfoFromEmbl")) setGeneric("addRepeatInfoFromEmbl", function(.Object, ...) standardGeneric("addRepeatInfoFromEmbl"))
#' addRepeatInfoFromEmbl-methods
#'
#' Add annotation contained in a RepBaseUpdate EMBL file to the reference of REs
#'
#' @param .Object         \code{\linkS4class{RepeatReference}} object
#' @param emblFile        Filename of the EMBL file from which the repeat annotaiton is retrieved. Can be NULL
#' @return modified \code{\linkS4class{RepeatReference}} object containing derived annotation
#'
#' @rdname addRepeatInfoFromEmbl-RepeatReference-method
#' @docType methods
#' @aliases addRepeatInfoFromEmbl
#' @aliases addRepeatInfoFromEmbl,RepeatReference-method
#' @author Fabian Mueller
#' @noRd
## @export
setMethod("addRepeatInfoFromEmbl", signature(.Object="RepeatReference"),
	function(.Object, emblFile=NULL){
		# check if there is an embl file corresponding to the reference fasta file
		if (is.null(emblFile)){
			require(tools)
			emblFile <- paste0(file_path_sans_ext(.Object@reference), ".embl")
		}
		if (file.exists(emblFile)){
			ri <- getRepeatFromEmbl(emblFile)
			repIds <- getRepeatIds(.Object)
			if (length(setdiff(repIds, names(ri)))>0){
				stop("Not all repeats are annotated in the embl file")
			}
			ri <- ri[repIds]
			.Object@repeatInfoList <- ri
		} else {
			logger.info(c("Could not find valid *.embl file at", emblFile))
		}
		return(.Object)
	}
)

if (!isGeneric("addRepeatInfoFromGenomeTrack")) setGeneric("addRepeatInfoFromGenomeTrack", function(.Object, ...) standardGeneric("addRepeatInfoFromGenomeTrack"))
#' addRepeatInfoFromGenomeTrack-methods
#'
#' Add annotation contained a genomic ranges track to the reference of REs
#'
#' @param .Object         \code{\linkS4class{RepeatReference}} object or \code{NULL} for default
#' @param grt    		  \code{\linkS4class{GenomeRepeatTrack}} object
#' @return modified \code{\linkS4class{RepeatReference}} object containing derived annotation
#'
#' @rdname addRepeatInfoFromGenomeTrack-RepeatReference-method
#' @docType methods
#' @aliases addRepeatInfoFromGenomeTrack
#' @aliases addRepeatInfoFromGenomeTrack,RepeatReference-method
#' @author Fabian Mueller
#' @noRd
## @export
setMethod("addRepeatInfoFromGenomeTrack", signature(.Object="RepeatReference"),
	function(.Object, grt=NULL){
		# check if there is an embl file corresponding to the reference fasta file
		if (is.null(grt)){
			genomeAss <- .config$assembly
			if (is.null(genomeAss)){
				if (.Object@species=="human"){
					genomeAss <- "hg38"
				} else if (.Object@species=="mouse"){
					genomeAss <- "mm10"
				} else {
					stop("Unknown species")
				}
			}
			grt <- GenomeRepeatTrack(genomeAss)
		}
		repIds <- getRepeatIds(.Object)
		numBases <- getRepeatGenomeCovg(grt)
		if (length(setdiff(repIds, names(numBases)))>0){
			stop("Not all repeats are annotated in the embl file")
		}		
		.Object@repeatInfo[repIds,"baseCovg"] <- numBases[repIds]
		genomeLength <- sum(seqlengths(getRepeatInstances(grt)))
		.Object@repeatInfo[repIds,"percCovg"] <- numBases[repIds]/genomeLength
		return(.Object)
	}
)


if (!isGeneric("filterRepeats_wl")) setGeneric("filterRepeats_wl", function(.Object, ...) standardGeneric("filterRepeats_wl"))
#' filterRepeats_wl-methods
#'
#' discard all but the specified REs from the repeat reference
#'
#' @param .Object         \code{\linkS4class{RepeatReference}} object
#' @param whitelist       identifiers or indices of REs to keep
#' @return modified \code{\linkS4class{RepeatReference}} object
#'
#' @rdname filterRepeats_wl-RepeatReference-method
#' @docType methods
#' @aliases filterRepeats_wl
#' @aliases filterRepeats_wl,RepeatReference-method
#' @author Fabian Mueller
#' @noRd
## @export
setMethod("filterRepeats_wl", signature(.Object="RepeatReference"),
	function(.Object, whitelist){
		res <- .Object
		isInWl <- .Object@repeatInfo$id %in% whitelist
		res@repeatInfo     <- res@repeatInfo[isInWl,]
		res@sequences      <- res@sequences[isInWl]
		res@repeatInfoList <- res@repeatInfoList[isInWl]
		return(res)
	}
)
# rr <- RepeatReference()
# keepRepIds <- sample(rr@repeatInfo$id, 100)
# rrf <- filterRepeats_wl(rr,keepRepIds)

if (!isGeneric("getSpecies")) setGeneric("getSpecies", function(.Object) standardGeneric("getSpecies"))
#' getSpecies-methods
#'
#' Retrieve the species the reference pertains to
#'
#' @param .Object         \code{\linkS4class{RepeatReference}} object
#' @return string specifying the species
#'
#' @rdname getSpecies-RepeatReference-method
#' @docType methods
#' @aliases getSpecies
#' @aliases getSpecies,RepeatReference-method
#' @author Fabian Mueller
#' @noRd
## @export
setMethod("getSpecies", signature(.Object="RepeatReference"),
	function(.Object){
		if (!.hasSlot(.Object, "species")){
			res <- .config$species
			logger.warning(c("The RepeatReference object is depreacted and does not have a species slot yet. --> using the corresponding config element (", res , ")"))
			return(res)
		}
		return(.Object@species)
	}
)
if (!isGeneric("getRepeatInfo")) setGeneric("getRepeatInfo", function(.Object) standardGeneric("getRepeatInfo"))
#' getRepeatInfo-methods
#'
#' Retrieve the repeat annotation data.frame from the repeat reference
#'
#' @param .Object         \code{\linkS4class{RepeatReference}} object
#' @return the \code{repeatInfo} as data.frame
#'
#' @rdname getRepeatInfo-RepeatReference-method
#' @docType methods
#' @aliases getRepeatInfo
#' @aliases getRepeatInfo,RepeatReference-method
#' @author Fabian Mueller
#' @noRd
## @export
setMethod("getRepeatInfo", signature(.Object="RepeatReference"),
	function(.Object){
		return(.Object@repeatInfo)
	}
)
if (!isGeneric("getRepeatIds")) setGeneric("getRepeatIds", function(.Object) standardGeneric("getRepeatIds"))
#' getRepeatIds-methods
#'
#' Retrieve the rpeat ids from the repeat reference
#'
#' @param .Object         \code{\linkS4class{RepeatReference}} object
#' @return character vector of repeat ids
#'
#' @rdname getRepeatIds-RepeatReference-method
#' @docType methods
#' @aliases getRepeatIds
#' @aliases getRepeatIds,RepeatReference-method
#' @author Fabian Mueller
#' @noRd
## @export
setMethod("getRepeatIds", signature(.Object="RepeatReference"),
	function(.Object){
		return(.Object@repeatInfo[,"id"])
	}
)
if (!isGeneric("getSequences")) setGeneric("getSequences", function(.Object) standardGeneric("getSequences"))
#' getSequences-methods
#'
#' Retrieve the repeat sequences from the reference of REs
#'
#' @param .Object         \code{\linkS4class{RepeatReference}} object
#' @return \code{DNAStringSet} containing sequence information for each RE.
#'
#' @rdname getSequences-RepeatReference-method
#' @docType methods
#' @aliases getSequences
#' @aliases getSequences,RepeatReference-method
#' @author Fabian Mueller
#' @noRd
## @export
setMethod("getSequences", signature(.Object="RepeatReference"),
	function(.Object){
		return(.Object@sequences)
	}
)
if (!isGeneric("getKmerCounts")) setGeneric("getKmerCounts", function(.Object, ...) standardGeneric("getKmerCounts"))
#' getKmerCounts-methods
#'
#' Count the number of all k-mers in the sequences of the REs
#'
#' @param .Object   \code{\linkS4class{RepeatReference}} object
#' @param k         length of the resulting k-mers
#' @return a matrix containing counts of all sequence k-mers (rows) for each RE (columns)
#'
#' @rdname getKmerCounts-RepeatReference-method
#' @docType methods
#' @aliases getKmerCounts
#' @aliases getKmerCounts,RepeatReference-method
#' @author Fabian Mueller
#' @noRd
## @export
setMethod("getKmerCounts", signature(.Object="RepeatReference"),
	function(.Object, k=4L){
		if (k<1L) stop("Invalid value for k")
		bases <- c("A","C","G","T")
		kmers <- c()
		for (i in 1:k){
			kmers.cur <- apply(expand.grid(rep(list(bases),i)), 1, FUN=function(x){paste(x,collapse="")})
			kmers <- c(kmers, kmers.cur)
		}
		kmerDSS <- DNAStringSet(kmers)

		kmerCounts <- do.call("cbind",lapply(.Object@sequences,FUN=function(ss){
			countPDict(kmerDSS, ss)
		}))
		colnames(kmerCounts) <- .Object@repeatInfo$id
		rownames(kmerCounts) <- kmers
		#dd <- dist(t(kmerCounts))
		return(kmerCounts)
	}
)

#' parseEMBL
#'
#' Parse an EMBL file
#'
#' @param fn file name
#' @return A list containing one element for each entry. Each entry in turn is a list containing a character vector of values for each field
#'
#' @author Fabian Mueller
#' @noRd
parseEMBL <- function(fn){
	# fn <- "~/tmp_work/repeats/epiRepeatR/repBaseDownload/Homo_sapiens_all.embl"
	require(stringr)
	elemSep <- "//"
	elemCont <- "  " #key for continuing current key
	txt <- readLines(fn)
	keys <- str_sub(txt, 1, 2)
	vals <- str_trim(str_sub(txt, 3))

	#assign keys for continuing fields: assign to previous key as long as there is a continuing key
	if (keys[1]==elemCont) stop("first line may not be continued fiels") #avoid endless loop if first line is continued field
	while(sum(keys==elemCont)){
		isContinued <- keys==elemCont
		isContinued[1] <- FALSE
		keys[isContinued] <- keys[which(isContinued)-1]
	}

	#remove empty value lines
	isEmptyLine <- vals=="" & keys != elemSep
	keys <- keys[!isEmptyLine]
	vals <- vals[!isEmptyLine]

	sepInds <- grep(elemSep, keys)
	#remove separator if its the first or last line in file
	if (is.element(1, sepInds)){
		sepInds <- sepInds[-1]
		keys <- keys[-1]
		vals <- vals[-1]
	}
	if (is.element(length(keys), sepInds)){
		sepInds <- sepInds[-length(keys)]
		keys <- keys[-length(keys)]
		vals <- vals[-length(keys)]
	}
	elemStarts <- c(1, sepInds+1)
	elemEnds   <- c(sepInds-1, length(keys))
	elems <- lapply(1:length(elemStarts), FUN=function(i){
		if (elemStarts[i] < elemEnds[i]){
			curKeys <- keys[elemStarts[i]:elemEnds[i]]
			curVals <- vals[elemStarts[i]:elemEnds[i]]
			curKeys.rle <- rle(curKeys)
			keyLengthsCum <- cumsum(curKeys.rle$lengths)
			keyStarts <- c(0, keyLengthsCum[-length(keyLengthsCum)]) + 1
			keyEnds <- keyLengthsCum
			res <- lapply(1:length(curKeys.rle$lengths), FUN=function(j){
				curVals[keyStarts[j]:keyEnds[j]]
			})
			names(res) <- curKeys.rle$values
			return(res)
		} else {
			return(NULL)
		}
	})
	isNullElem <- vapply(elems, is.null, logical(1))
	elems <- elems[!isNullElem]
	return(elems)
}

#' getRepeatFromEmbl
#'
#' Parse an EMBL file and return a list of repeat elements along with their annotation
#'
#' @param fn file name (EMBL file)
#' @return A list containing one element for each repeat element. Each entry in turn is a list containing a character vector of values for each annotation
#'
#' @author Fabian Mueller
#' @noRd
getRepeatFromEmbl <- function(fn){
	require(stringr)
	getTerms <- function(x){
		r <- strsplit(paste(x,collapse=";"), ";")[[1]]
		r <- str_trim(r)
		#remove leading and trailing punctation characters
		r <- sub("^\\.+","",r)
		r <- sub("\\.+$","",r)
		r <- r[nchar(r)>0]
		return(r)
	}
	getSeq <- function(x){
		if (grepl("^Sequence",x[1])) x <- x[-1]
		x <- paste(x, collapse="")
		#remove everything that is not base nomenclature
		x <- gsub("[^acgtnACGTNrywsmkbhdvRYWSMKBHDV]", "", x)
	}
	elems <- parseEMBL(fn)
	repList <- lapply(elems, FUN=function(x){
		idLine <- strsplit(x[["ID"]][1],";")[[1]]
		curName <- strsplit(idLine[1]," ")[[1]][1]
		curLvlAbr <- str_trim(idLine[3]) #species level abreviation from ID field
		curSpecies <- x[["OS"]][1]
		curSpeciesTerms <- getTerms(x[["OC"]])
		curRepeatTerms  <- getTerms(x[["KW"]])
		curRepeatFamily <- curRepeatTerms[1]
		curSeq <- DNAString(getSeq(x[["SQ"]]))
		res <- list(
			name=curName,
			family=curRepeatFamily,
			repeatTerms=curRepeatTerms,
			species=curSpecies,
			speciesTerms=curSpeciesTerms,
			specLvlAbr=curLvlAbr,
			sequence=curSeq,
			asParsed=x
		)
		return(res)
	})
	repeatIds <- vapply(repList,FUN=function(x){x$name},character(1))
	if (any(duplicated(repeatIds))) stop("EMBL file contains duplicated repeat ids")
	names(repList) <- repeatIds
	# class(repList) <- c("RepeatAnnotation", "RepeatAnnotationEmbl")
	return(repList)
}
