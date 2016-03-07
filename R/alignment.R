#' indexRepeatReference.bwa
#'
#' Index repeat reference file for 'bwa'
#'
#' @param refFn			repeat reference to align to (fasta)
#' @return invisible(0) if successfull
#'
#' @author Fabian Mueller
#' @noRd
indexRepeatReference.bwa <- function(refFn=.config$refFasta){
	res <- 0
	indexExt <- c("amb","ann","bwt","pac","sa")# file extensions required for a reference file to be indexed
	isIndexed <- all(file.exists(paste0(refFn,".",indexExt)))
	if (!isIndexed){
		cmd <- paste(
			"bwa index",
			refFn
		)
		logger.info(c("Indexing reference for BWA alignment:",refFn))
		errCode <- system(cmd)
		if (errCode != 0) logger.error(c("Command line call failed. Command used:", cmd))
	}
	invisible(res)
}
