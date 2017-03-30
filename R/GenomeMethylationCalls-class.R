setClassUnion("ListOrNULL", c("list", "NULL"))

setClass("GenomeMethylationCalls",
	slots = list(
		repTrack   = "GenomeRepeatTrack",
		methCalls  = "ListOrNULL"
	),
	package = "epiRepeatR"
)

setMethod("initialize", "GenomeMethylationCalls",
	function(.Object,
		mcGr,
		repeatTrack
	){
		repTrack <- getRepeatInstances(repeatTrack)

		seqlvls.stripped <- gsub("^chr", "", seqlevels(repTrack))

		if (sum(seqlevels(mcGr) %in% seqlvls.stripped) > sum(seqlevels(mcGr) %in% seqlevels(repTrack))){
			logger.info("Stripping 'chr' prefix from reference repeat coordinates")
			seqlevels(repTrack) <- seqlvls.stripped
		}

		seqlvls.both <- intersect(seqlevels(repTrack), seqlevels(mcGr))
		if (length(seqlvls.both) < length(seqlevels(repTrack))){
			excl.seqlvls <- setdiff(seqlevels(repTrack), seqlvls.both)
			logger.warning(c("The following chroms are not present in the methylation calls and are disregarded in the analysis:", paste(excl.seqlvls, collapse=",")))
			repTrack <- endoapply(repTrack, FUN=function(x){
				x[!(as.character(seqnames(x)) %in% excl.seqlvls)]
			})
			seqlevels(repTrack) <- setdiff(seqlevels(repTrack), excl.seqlvls)
		}
		if (length(seqlvls.both) < length(seqlevels(mcGr))){
			excl.seqlvls <- setdiff(seqlevels(mcGr), seqlvls.both)
			logger.warning(c("The following chroms are not present in the repeat annotation and are disregarded in the analysis:", paste(excl.seqlvls, collapse=",")))
			mcGr <- mcGr[!(as.character(seqnames(mcGr)) %in% excl.seqlvls)]
			seqlevels(mcGr) <- setdiff(seqlevels(mcGr), excl.seqlvls)
		}

		logger.start("Computing overlaps for repeats")
			mcs <- lapply(repTrack, FUN=function(repGr){
				oo <- as.matrix(findOverlaps(mcGr, repGr, ignore.strand=TRUE))
				mcInd <- unique(oo[,1])
				repInds <- tapply(oo[,2], oo[,1], FUN=c, simplify=FALSE)
				names(repInds) <- NULL
				mcGr.sub <- mcGr[mcInd]
				emd.sub <- elementMetadata(mcGr.sub)
				methCalls <- data.frame(
					cPos=start(mcGr.sub),
					numM=emd.sub$M,
					numU=emd.sub$T-emd.sub$M,
					numT=emd.sub$T,
					repeatInstanceIndex=repInds,
					cChrom=as.character(seqnames(mcGr.sub)),
					cStrand=as.character(strand(mcGr.sub))
				)
				res <- list(
					methCalls=methCalls,
					readStats=c("numReads"=NA, "numReads_informative"=NA)
				)
				return(res)
			})
		logger.completed()

		.Object@repTrack  <- repeatTrack
		.Object@methCalls <- mcs
		.Object
	}
)
GenomeMethylationCalls <- function(mcFile, repeatTrack, mcFormat=getConfigElement("meth.methCallFormat")){
	logger.start("Reading methylation calls")
	if (mcFormat=="BisSNP"){
		mcGr <- parseMcTable.bissnp(mcFile)
	} else if (mcFormat=="EPP"){
		mcGr <- parseMcTable.epp(mcFile)
	} else {
		logger.error(c("Unrecognized methylation call format:", mcFormat))
	}
	logger.completed()

	object <- new("GenomeMethylationCalls",
		mcGr, repeatTrack
	)
	return(object)
}

if (!isGeneric("getMethylationCalls")) setGeneric("getMethylationCalls", function(.obj, ...) standardGeneric("getMethylationCalls"))
setMethod("getMethylationCalls", signature(.obj="GenomeMethylationCalls"),
	function(.obj, ...){
		mcs <- .obj@methCalls
		class(mcs) <- c("RepeatMethCalls")
		return(mcs)
	}
)
