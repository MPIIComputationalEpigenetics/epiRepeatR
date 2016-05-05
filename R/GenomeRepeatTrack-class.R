setClass("GenomeRepeatTrack",
	slots = list(
		genome = "character",
		repInstances = "GRangesList"
	),
	package = "epiRepeatR"
)

setMethod("initialize", "GenomeRepeatTrack",
	function(.Object,
		genome,
		trackFile=NULL
	){
		columnMatching <- c(
			chrom = "genoName",
			start = "genoStart",
			end = "genoEnd",
			strand = "strand",
			repName = "repName",
			repClass = "repClass",
			repFamily = "repFamily"

		)
		if (is.null(trackFile)){
			require(rtracklayer)
			session <- browserSession()
			genome(session) <- genome
			tt <- getTable(ucscTableQuery(session, track="RepeatMasker", table="rmsk"))
		} else {
			tt <- read.table(trackFile, sep="\t", header=TRUE, stringsAsFactors=FALSE, comment.char="")
		}
		missingCols <- !(columnMatching %in% colnames(tt))
		if (any(missingCols)){
			stop(paste0("The following columns were not found in the repeat track annotation", paste(columnMatching[missingCols], collapse=",")))
		}
		repeatTrackGr <- GRanges(
			seqnames=tt[,columnMatching["chrom"]],
			ranges=IRanges(start=tt[,columnMatching["start"]], end=tt[,columnMatching["end"]]),
			strand=tt[,columnMatching["strand"]],
			seqinfo=Seqinfo(genome=genome)
		)
		elementMetadata(repeatTrackGr) <- tt[,columnMatching[c("repName", "repClass", "repFamily")]]
		reps <- elementMetadata(repeatTrackGr)[,"repName"]
		repNames <- sort(unique(reps))
		repInstances <- lapply(repNames, FUN=function(rn){
			repeatTrackGr[reps==rn]
		})
		names(repInstances) <- repNames
		repInstances <- GRangesList(repInstances)

		.Object@genome <- genome
		.Object@repInstances <- repInstances
		.Object
	}
)
GenomeRepeatTrack <- function(genome, trackFile=NULL){
	object <- new("GenomeRepeatTrack",
		genome, trackFile
	)
	return(object)
}

if (!isGeneric("getRepeatInstances")) setGeneric("getRepeatInstances", function(.obj) standardGeneric("getRepeatInstances"))
setMethod("getRepeatInstances", signature(.obj="GenomeRepeatTrack"),
	function(.obj){
		return(.obj@repInstances)
	}
)
