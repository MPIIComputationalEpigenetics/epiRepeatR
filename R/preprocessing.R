#' getScanBamFilterFlags
#'
#' Convenience function to extract filtering flags for the \code{extractFastqFromBam} function.
#'
#' @param onlyQCpass			Only reads which pass the QC filter should be returned
#' @param mappingStatus			One of 'all', 'mapped' (returns only mapped reads), 'unmapped' (returns only unmapped reads)
#' @param pairedReads			One of 'all', 'read1' (returns only the first read of each pair), 'read2' (returns only the second read of each pair)
#' @return a filtering object as returned by \code{Rsamtools::scanBamFlag}
#'
#' @author Fabian Mueller
#' @noRd
getScanBamFilterFlags <- function(onlyQCpass=TRUE, mappingStatus="all", pairedReads="all"){
	bam.flag.isPaired <- NA
	bam.flag.isUnmappedQuery <- NA
	bam.flag.isFirstMateRead <- NA
	bam.flag.isSecondMateRead <- NA
	bam.flag.isNotPassingQualityControls <- NA
	if (onlyQCpass){
		bam.flag.isNotPassingQualityControls <- FALSE
	}
	if (mappingStatus=="all"){
		bam.flag.isUnmappedQuery  <- NA
	} else if (mappingStatus=="mapped"){
		bam.flag.isUnmappedQuery <- FALSE
	} else if (mappingStatus=="unmapped"){
		bam.flag.isUnmappedQuery <- TRUE
	} else {
		logger.error("invalid value for mappingStatus parameter. must be one of 'all','mapped','unmapped'.")
	}
	if (pairedReads=="all"){
		bam.flag.isFirstMateRead  <- NA
		bam.flag.isSecondMateRead <- NA
	} else if (pairedReads=="read1"){
		bam.flag.isFirstMateRead  <- TRUE
		bam.flag.isSecondMateRead <- FALSE
	} else if (pairedReads=="read2"){
		bam.flag.isFirstMateRead  <- FALSE
		bam.flag.isSecondMateRead <- TRUE
	} else {
		logger.error("invalid value for 'pairedReads' parameter. must be one of 'all','read1','read2'.")
	}
	res <- scanBamFlag(
		isPaired=bam.flag.isPaired,
		isUnmappedQuery=bam.flag.isUnmappedQuery,
		isFirstMateRead=bam.flag.isFirstMateRead,
		isSecondMateRead=bam.flag.isSecondMateRead,
		isNotPassingQualityControls=bam.flag.isNotPassingQualityControls
	)
	return(res)
}
# #' extractBamFromBam
# #'
# #' Extract reads from bam files and save them to bam
# #'
# #' @param inputBam				input bam file name
# #' @param outputBam				output bam file name
# #' @param scanBamFilterFlags	read flags to be filtered for during the extraction process. as returned by \code{getScanBamFilterFlags} or \code{Rsamtools::scanBamFlag}
# #' @return a list with details on the extraction process (e.g. the status of success of failure)
# #'
# #' @author Fabian Mueller
# #' @noRd
# extractBamFromBam <- function(inputBam, outputBam, scanBamFilterFlags=getScanBamFilterFlags()){
# 	filterBam(inputBam, outputBam, param=ScanBamParam(flag=scanBamFilterFlags))
# 	# leads to memory errors:
# 	#  *** caught segfault ***
# 	# address (nil), cause 'memory not mapped'
# 	# --> maybe resort to command line samtools

# 	res <- list(status="success")
# 	return(res)
# }
