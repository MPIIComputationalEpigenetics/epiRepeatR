#' runAnalysis
#'
#' Run the epiRepeatR pipeline. This is the main function.
#'
#' @param anaDir	  Analysis directory. epiRepeatR will manage this directory. If it is non-existing a new analysis will be launched
#'					  if it is existing, epiRepeatR will try to resume the analysis
#' @param fileTable	  Table containing file paths of the input files and annotation data. Will be ignored, when resuming an analysis.
#' @param resetToStep If resuming an existing analysis or running from an existing analysis directory, reset the analysis
#'                    to the beginning of the specified analysis step prior to running the pipeline
#' @param submission  Type of submission to be used. "system" (default) for system calls "sge" for Sun Grid Engine.
#' @param cmdrArgs    Named list of additional arguments for the construction of the CommandR object used for executing the pipeline. 
#' @return nothing of particular interest
#'
#' @author Fabian Mueller
#' @export
runAnalysis <- function(anaDir, fileTable=NULL, resetToStep=NULL, submission="system", cmdrArgs=list()){
	logger.start(fname=NA)
	newAna <- TRUE
	logger.start("Analysis preliminaries")
		logger.info(c("Output directory:",anaDir))
		piprObjPath <- file.path(anaDir, "status", "pipr.rds")
		if (dir.exists(anaDir)) {
			newAna <- FALSE
			logger.info("Output directory already exists --> continuing analysis")
			if (!file.exists(piprObjPath)){
				logger.error("Analysis directory exists but does not contain pipeline object")
			}
			if (!is.null(fileTable)){
				logger.warning("fileTable argument suplied but not required")
			}
			logger.status("Loading pipeline object from directory...")
			pipr <- readRDS(piprObjPath)
		} else {
			if (is.null(fileTable)){
				logger.error("fileTable argument required")
			}
			logger.status("Reading files and annotation...")
			anaMan <- readFileTable(fileTable)
			logger.start("Constructing analysis pipeline")
				pipr <- buildPipeline(anaMan, anaDir)
			logger.completed()
		}
		usePiprTempDir <- TRUE
		if (usePiprTempDir) .config$tempDir <-  getDir(pipr,"temp")
	logger.completed()
	logger.start(fname=c(NA,file.path(getDir(pipr,"log"), format(Sys.time(), "run_%Y%m%d_%H_%M_%S.log"))))

	if (!newAna && !is.null(resetToStep)){
		pipr <- resetStep(pipr, resetToStep)
		# overwrite existing configuration file
		cfgDir <- file.path(getDir(pipr, "base"), "config")
		cfgSavePath <- file.path(cfgDir, "config.json")
		if (file.exists(cfgSavePath)){
			logger.warning(c("Overwriting saved configuration file @",cfgSavePath))
		}
		saveConfig(cfgSavePath)
	}
	if (submission == "system"){
		cmdr <- do.call("CommandRsystem", c(list(logDir=getDir(pipr, "log")), cmdrArgs))
	} else if (submission == "sge"){
		cmdr <- do.call("CommandRsge", c(list(logDir=getDir(pipr, "log")), cmdrArgs))
	} else if (submission == "slurm"){
		cmdr <- do.call("CommandRslurm", c(list(logDir=getDir(pipr, "log")), cmdrArgs))
	} else {
		logger.error("Could not execute pipeline: unknown submission type.")
	}

	run(pipr, cmdr, logCommands=TRUE)
}
