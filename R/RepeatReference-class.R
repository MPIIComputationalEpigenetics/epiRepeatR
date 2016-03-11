setClassUnion("ListOrNULL", c("list", "NULL"))

setClass("RepeatReference",
	slots = list(
		reference="character",
		repeatInfo="data.frame",
		repeatInfoList="ListOrNULL",
		sequences="DNAStringSet"
	),
	package = "epiRepeatR"
)

setMethod("initialize", "RepeatReference",
	function(.Object,
		reference=.config$refFasta
	){
		.Object@reference=reference
		.Object@repeatInfo=getReferenceInfo(reference)
		.Object@repeatInfoList=NULL
		.Object@sequences=getDNAStringForReference(reference)

		.Object
	}
)

RepeatReference <- function(reference=.config$refFasta){
	obj <- new("RepeatReference",
		reference
	)
	obj <- addRepeatInfoFromEmbl(obj)
	return(obj)
}


getDNAStringForReference <- function(reference){
	refSeqs <- readBStringSet(reference)
	#replace letters "x"
	for (i in 1:length(refSeqs)){
		refSeqs[[i]][start(matchPattern("x",refSeqs[[i]]))] <- "n"
	}
	refSeqs <- DNAStringSet(refSeqs)
	return(refSeqs)
}

getReferenceInfo <- function(reference){
	refSeqs <- readBStringSet(reference)
	seqLens <- vapply(refSeqs,length,integer(1))
	res <- data.frame(t(data.frame(strsplit(names(refSeqs),"\t"))),seqLength=seqLens)
	colnames(res)[1:4] <- c("id","family","species","seqLength")
	rownames(res) <- res[,1]
	res[,"id"] <- as.character(res[,"id"])
	return(res)
}

if (!isGeneric("addRepeatInfoFromEmbl")) setGeneric("addRepeatInfoFromEmbl", function(.Object, ...) standardGeneric("addRepeatInfoFromEmbl"))
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

if (!isGeneric("filterRepeats_wl")) setGeneric("filterRepeats_wl", function(.Object, ...) standardGeneric("filterRepeats_wl"))
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
if (!isGeneric("getRepeatInfo")) setGeneric("getRepeatInfo", function(.Object) standardGeneric("getRepeatInfo"))
setMethod("getRepeatInfo", signature(.Object="RepeatReference"),
	function(.Object){
		return(.Object@repeatInfo)
	}
)
if (!isGeneric("getRepeatIds")) setGeneric("getRepeatIds", function(.Object) standardGeneric("getRepeatIds"))
setMethod("getRepeatIds", signature(.Object="RepeatReference"),
	function(.Object){
		return(.Object@repeatInfo[,"id"])
	}
)
if (!isGeneric("getSequences")) setGeneric("getSequences", function(.Object) standardGeneric("getSequences"))
setMethod("getSequences", signature(.Object="RepeatReference"),
	function(.Object){
		return(.Object@sequences)
	}
)
if (!isGeneric("getKmerCounts")) setGeneric("getKmerCounts", function(.Object, ...) standardGeneric("getKmerCounts"))
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
