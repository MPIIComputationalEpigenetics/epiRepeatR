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

#' indexRepeatReference.bowtie2
#'
#' Index repeat reference file for 'bowtie2'
#'
#' @param refFn			repeat reference to align to (fasta)
#' @return invisible(0) if successfull
#'
#' @author Fabian Mueller
#' @noRd
indexRepeatReference.bowtie2 <- function(refFn=.config$refFasta){
	res <- 0
	indexExt <- c("1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2")# file extensions required for a reference file to be indexed
	isIndexed <- all(file.exists(paste0(refFn,".",indexExt)))
	if (!isIndexed){
		wdSave <- getwd()
		setwd(dirname(refFn))
		cmd <- paste(
			"bowtie2-build",
			basename(refFn),
			basename(refFn)
		)
		logger.info(c("Indexing reference for bowtie2 alignment:",refFn))
		errCode <- system(cmd)
		setwd(wdSave)
		if (errCode != 0) logger.error(c("Command line call failed. Command used:", cmd))
	}
	invisible(res)
}
