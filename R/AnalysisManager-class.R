#' AnalysisManager Class
#'
#' A class for managing a repeat epigenome analysis
#' 
#' @details
#' Provides everything needed to keep track of Sample annotation, generated files and steps to be taken in the analysis
#'
#' @section Slots:
#' \describe{
#'   \item{\code{sampleAnnot}}{Sample annotation table for the analysis}
#'   \item{\code{fileTable}}{table of all input, intermediate and output files of the analysis and their links to analysis steps} 
#'   \item{\code{stepDetails}}{details on steps carried out during the analysis}
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{\code{\link{getSampleAnnot,AnalysisManager-method}}}{retrieve the sample annotation information}
#'   \item{\code{\link{getFileTable,AnalysisManager-method}}}{retrieve the sample annotation information}
#'   \item{\code{\link{getSampleNames,AnalysisManager-method}}}{retrieve the sample annotation information}
#'   \item{\code{\link{getDataTypes,AnalysisManager-method}}}{retrieve the sample annotation information}
#'   \item{\code{\link{getMarks,AnalysisManager-method}}}{retrieve the sample annotation information}
#'   \item{\code{\link{getSampleMarkDatatypes,AnalysisManager-method}}}{retrieve the sample annotation information}
#'   \item{\code{\link{addFile,AnalysisManager-method}}}{retrieve the sample annotation information}
#'   \item{\code{\link{getStructure4analysisStep,AnalysisManager-method}}}{retrieve the sample annotation information}
#'   \item{\code{\link{setStepDetails,AnalysisManager-method}}}{retrieve the sample annotation information}
#'   \item{\code{\link{prependAnalysisDirForFilename,AnalysisManager-method}}}{retrieve the sample annotation information}
#' }
#'
#' @noRd
#' @name RnBDiffMeth-class
#' @rdname RnBDiffMeth-class
#' @author Fabian Mueller
setClass("AnalysisManager",
	slots = list(
		sampleAnnot = "data.frame",
		fileTable = "data.frame"
	),
	package = "epiRepeatR"
)

#' initialize.AnalysisManager
#'
#' Initialize an AnalysisManager object
#' 
#' @param .Object New instance of \code{AnalysisManager}.
#' @param fName Filename of the Analysis configuration table. It should contain columns for filenames of the input files (default: fileName), the data type of the input files (default column name: dataType), a sample identifier (SampleName), an epigenetic mark identifier column (mark) and a column indicating the analysis Step that generated the input file (analysisStep)
#' @param groupColumns  column names of sample annotation
#' @param columnMap   mapping of column names to column types (see details on \code{fName} for default column names)
#'
#' @noRd
#' @author Fabian Mueller
#' @docType methods
setMethod("initialize", "AnalysisManager",
	function(.Object,
		fName,
		groupColumns=c(),
		columnMap=c(
			fileName="fileName",
			dataType="dataType",
			sampleName="sampleName",
			mark="mark",
			analysisStep="analysisStep"
		)
	){
		ft <- read.table(fName, sep="\t", header=TRUE, stringsAsFactors=FALSE)
		ft.base <- data.frame(ft[,columnMap[c("fileName","dataType","mark","sampleName","analysisStep")]], stringsAsFactors=FALSE)
		colnames(ft.base) <- c("fileName","dataType","mark","sampleName","analysisStep")

		ft.base <- data.frame(ft.base, fileType=getFileTypeFromFileName(ft.base$fileName), stringsAsFactors=FALSE)

		sampleNames <- sort(unique(ft.base[,"sampleName"]))
		
		sampleAnnot <- data.frame(matrix(nrow=length(sampleNames),ncol=0))

		for (grp in groupColumns){
			grpAnnot <- ft[,grp]
			curAnnot <- sapply(sampleNames,FUN=function(sn){
				vals <- grpAnnot[ft.base[,"sampleName"]==sn]
				if (length(unique(vals))!=1){
					logger.warning(c("multiple non-matching annotations found for sample",sn,"and annotation",grp,"--> picking the first one"))
				}
				return(vals[1])
			})
			sampleAnnot[[grp]] <- curAnnot
		}
		rownames(sampleAnnot) <- sampleNames

		.Object@sampleAnnot <- sampleAnnot
		.Object@fileTable <- ft.base
		.Object
	}
)

#' getSampleAnnot-methods
#'
#' Retrieve the sample annotation table from a \code{\linkS4class{AnalysisManager}} object
#'
#' @param .Object \code{\linkS4class{AnalysisManager}} object
#' @return The sample annotation table from the object
#'
#' @rdname getSampleAnnot-AnalysisManager-methods
#' @docType methods
#' @author Fabian Mueller
#' @aliases getSampleAnnot,AnalysisManager-method
#' @noRd
if (!isGeneric("getSampleAnnot")) setGeneric("getSampleAnnot", function(.Object) standardGeneric("getSampleAnnot"))
setMethod("getSampleAnnot", signature(.Object="AnalysisManager"),
	function(.Object){
		return(.Object@sampleAnnot)
	}
)
#' getFileTable-methods
#'
#' Retrieve the file annotation table from a \code{\linkS4class{AnalysisManager}} object
#'
#' @param .Object \code{\linkS4class{AnalysisManager}} object
#' @return The file annotation table from the object
#'
#' @rdname getFileTable-AnalysisManager-methods
#' @docType methods
#' @author Fabian Mueller
#' @aliases getFileTable,AnalysisManager-method
#' @noRd
if (!isGeneric("getFileTable")) setGeneric("getFileTable", function(.Object) standardGeneric("getFileTable"))
setMethod("getFileTable", signature(.Object="AnalysisManager"),
	function(.Object){
		return(.Object@fileTable)
	}
)
#' getSampleNames-methods
#'
#' Retrieve a list of all sample names for the analysis described in an \code{\linkS4class{AnalysisManager}} object
#'
#' @param .Object \code{\linkS4class{AnalysisManager}} object
#' @return a vector of sample names (sorted alphanumerically)
#'
#' @rdname getSampleNames-AnalysisManager-methods
#' @docType methods
#' @author Fabian Mueller
#' @aliases getSampleNames,AnalysisManager-method
#' @noRd
if (!isGeneric("getSampleNames")) setGeneric("getSampleNames", function(.Object) standardGeneric("getSampleNames"))
setMethod("getSampleNames", signature(.Object="AnalysisManager"),
	function(.Object){
		return(sort(unique(.Object@fileTable[,"sampleName"])))
	}
)
#' getDataTypes-methods
#'
#' Retrieve a list of all data types for the analysis described in an \code{\linkS4class{AnalysisManager}} object
#'
#' @param .Object \code{\linkS4class{AnalysisManager}} object
#' @return a vector of data types (sorted alphanumerically)
#'
#' @rdname getDataTypes-AnalysisManager-methods
#' @docType methods
#' @author Fabian Mueller
#' @aliases getDataTypes,AnalysisManager-method
#' @noRd
if (!isGeneric("getDataTypes")) setGeneric("getDataTypes", function(.Object) standardGeneric("getDataTypes"))
setMethod("getDataTypes", signature(.Object="AnalysisManager"),
	function(.Object){
		return(sort(unique(.Object@fileTable[,"dataType"])))
	}
)
#' getMarks-methods
#'
#' Retrieve a list of all epigenetic marks for the analysis described in an \code{\linkS4class{AnalysisManager}} object
#'
#' @param .Object \code{\linkS4class{AnalysisManager}} object
#' @return a vector of epigenetic marks (sorted alphanumerically)
#'
#' @rdname getMarks-AnalysisManager-methods
#' @docType methods
#' @author Fabian Mueller
#' @aliases getMarks,AnalysisManager-method
#' @noRd
if (!isGeneric("getMarks")) setGeneric("getMarks", function(.Object) standardGeneric("getMarks"))
setMethod("getMarks", signature(.Object="AnalysisManager"),
	function(.Object){
		return(sort(unique(.Object@fileTable[,"mark"])))
	}
)
#' getSampleMarks-methods
#'
#' Retrieve a list of epigenetic marks for a given sample id in the analysis described in an \code{\linkS4class{AnalysisManager}} object
#'
#' @param .Object \code{\linkS4class{AnalysisManager}} object
#' @param sampleName sample id
#' @return a vector of epigenetic marks (sorted alphanumerically)
#'
#' @rdname getSampleMarks-AnalysisManager-methods
#' @docType methods
#' @author Fabian Mueller
#' @aliases getSampleMarks,AnalysisManager-method
#' @noRd
if (!isGeneric("getSampleMarks")) setGeneric("getSampleMarks", function(.Object, ...) standardGeneric("getSampleMarks"))
setMethod("getSampleMarks", signature(.Object="AnalysisManager"),
	function(.Object, sampleName){
		return(sort(unique(.Object@fileTable[.Object@fileTable[,"sampleName"]==sampleName,"mark"])))
	}
)
#' getSampleMarkDatatypes-methods
#'
#' Retrieve a list of all data types for a given sample id and epigenetic mark in the analysis described in an \code{\linkS4class{AnalysisManager}} object
#'
#' @param .Object \code{\linkS4class{AnalysisManager}} object
#' @param sampleName sample id
#' @param mark epigenetic mark id
#' @return a vector of data types (sorted alphanumerically)
#'
#' @rdname getSampleMarkDatatypes-AnalysisManager-methods
#' @docType methods
#' @author Fabian Mueller
#' @aliases getSampleMarkDatatypes,AnalysisManager-method
#' @noRd
if (!isGeneric("getSampleMarkDatatypes")) setGeneric("getSampleMarkDatatypes", function(.Object, ...) standardGeneric("getSampleMarkDatatypes"))
setMethod("getSampleMarkDatatypes", signature(.Object="AnalysisManager"),
	function(.Object, sampleName, mark){
		return(sort(unique(.Object@fileTable[.Object@fileTable[,"sampleName"]==sampleName & .Object@fileTable[,"mark"]==mark,"dataType"])))
	}
)

#' addFile
#'
#' add a file to the file annotation table
#'
#' @param fileTable table of files as returned by \code{getFileTable}
#' @param fileName filename associated with the file
#' @param dataType data type associated with the file
#' @param mark epigenetic mark associated with the file
#' @param sampleName sample id associated with the file
#' @param analysisStep analysis step producing the file
#' @return the modified file table
#'
#' @author Fabian Mueller
#' @noRd
addFile <- 	function(fileTable, fileName, dataType, mark, sampleName, analysisStep){
	addRow <- data.frame(
		fileName=fileName,
		dataType=dataType,
		mark=mark,
		sampleName=sampleName,
		fileType=getFileTypeFromFileName(fileName),
		analysisStep=analysisStep,
		stringsAsFactors=FALSE
	)
	fileTable <- rbind(fileTable, addRow)
	return(fileTable)
}

#' buildPipeline-methods
#'
#' Construct a pipeline from the analysis description using the \code{muPipeR} package
#'
#' @param .Object \code{\linkS4class{AnalysisManager}} object
#' @return \code{\linkS4class{PipR}} object
#'
#' @rdname buildPipeline-AnalysisManager-methods
#' @docType methods
#' @author Fabian Mueller
#' @aliases buildPipeline,AnalysisManager-method
#' @noRd
if (!isGeneric("buildPipeline")) setGeneric("buildPipeline", function(.Object, ...) standardGeneric("buildPipeline"))
setMethod("buildPipeline", signature(.Object="AnalysisManager"),
	function(.Object, baseDir){
		# TODO:
		# - avoid full paths to scripts in inst/exdata --> make execution environment independent
		analysisSteps <- c(
			"seqReads", "bamExtract", "repeatAlignment",
			"methCalling", "chipQuantification",
			"plotRepeatGroupTreesMeth", "plotRepeatMarkTree")


		ft <- getFileTable(.Object)
		sampleNames <- getSampleNames(.Object)
		markNames <- getMarks(.Object)
		dataTypes <- getDataTypes(.Object)
		sampleMarks <- lapply(sampleNames, FUN=function(sn){
			getSampleMarks(.Object, sn)
		})
		names(sampleMarks) <- sampleNames

		# check if multiple analysis levels for the same combination are present
		for (sn in sampleNames){
			for (mn in sampleMarks[[sn]]){
				for (dn in getSampleMarkDatatypes(.Object, sn, mn)){
					expInds <- ft[,"sampleName"]==sn & ft[,"mark"]==mn & ft[,"dataType"]==dn
					if (length(unique(ft[expInds,"analysisStep"]))>1){
						logger.error(c("Multiple analysis steps found for", sn, mn, dn))
					}
				}
			}
		}

		#initialize pipeline object
		pipr <- PipR(baseDir)
		cfgDir         <- file.path(getDir(pipr, "base"), "config")
		cfgDir.forPipe <- file.path("${BASEDIR}","config")
		dir.create(cfgDir)

		cfgSavePath         <- file.path(cfgDir, "config.json")
		cfgSavePath.forPipe <- file.path(cfgDir.forPipe, "config.json")
		saveConfig(cfgSavePath)
		
		logger.info(c("Saved config file to",cfgSavePath))
		anamanSavePath <- file.path(cfgDir, "anaMan.rds")
		anamanSavePath.forPipe <- file.path(cfgDir.forPipe, "anaMan.rds")
		saveRDS(.Object, anamanSavePath)
		logger.info(c("Saved analysis manager to",anamanSavePath))
		fileTabPath <- file.path(cfgDir, "fileTable.tsv")

		tempDir <- getDir(pipr, "temp")

		# add analysis steps
		#-----------------------------------------------------------------------
		stepName <- "bamExtract"
		logger.status(c("step:", stepName))
		cmd.bamExtract <- "bash"
		bamExtractScript <- system.file(file.path("extdata", "exec", "bamExtract.sh"), package="epiRepeatR")
		args.bamExtract <- list()
		for (sn in sampleNames){
			for (mn in sampleMarks[[sn]]){
				for (dn in getSampleMarkDatatypes(.Object, sn, mn)){
					stepId <- paste(sn,mn,dn,stepName,sep="_")
					expInds <- ft[,"sampleName"]==sn & ft[,"mark"]==mn & ft[,"dataType"]==dn
					stepInds.input <- expInds & ft[,"fileType"]=="bam" & ft[,"analysisStep"]=="seqReads"
					inputPresent <- any(stepInds.input)
					doStep <- inputPresent
					if (doStep){
						if (sum(stepInds.input)>1){
							logger.warning(c("Multiple input bam files found for step",stepId,"--> picking the first one"))
						}
						inBamFn  <- ft[which(stepInds.input)[1],"fileName"]
						outBamFn <- file.path(paste0("${STEPDIR:", stepName, "}"), paste0(stepId,".bam"))

						#flags as returned by scanBamFlag (Rsamtools)
						sbf <- getScanBamFilterFlags(mappingStatus=.config$inputBam.mappingStatus)
						#need to be converted to fit with the "-f" and "-F" parameters of samtools
						arg.f <- bitwCompl(sbf["keep0"])
						arg.F <- bitwCompl(sbf["keep1"])

						curArgs <- unname(c(
							bamExtractScript,
							.config$samtools.exec,
							inBamFn,
							outBamFn,
							arg.f,
							arg.F
						))

						args.bamExtract <- c(args.bamExtract, list(curArgs))

						ft <- addFile(ft, outBamFn, dn, mn, sn, stepName)
					}
				}
			}
		}
		if (length(cmd.bamExtract) > 0 && length(args.bamExtract) > 0){
			pipr <- addStep(pipr, stepName, cmd.bamExtract, args.bamExtract, parents=character())
		}
		#-----------------------------------------------------------------------
		stepName <- "repeatAlignment"
		logger.status(c("step:", stepName))
		cmd.repeatAlignment <- "bash"
		args.repeatAlignment <- list()
		do.bwaIndex <- FALSE
		for (sn in sampleNames){
			for (mn in sampleMarks[[sn]]){
				for (dn in getSampleMarkDatatypes(.Object, sn, mn)){
					stepId <- paste(sn,mn,dn,stepName,sep="_")
					expInds <- ft[,"sampleName"]==sn & ft[,"mark"]==mn & ft[,"dataType"]==dn
					stepInds.input <- expInds & ft[,"fileType"]=="bam" & ft[,"analysisStep"]=="bamExtract"
					inputPresent <- any(stepInds.input)
					alnType <- NULL
					if (is.element(dn, c("WGBS","RRBS"))){
						alnType <- .config$aligner.bs
					} else if (is.element(dn,c("Input","ChIPseq"))){
						alnType <- .config$aligner.chip
					}
					curScript <- NULL
					curExec <- NULL
					curOtherArgs <- c()
					if (alnType == "bsmap"){
						curScript <- system.file(file.path("extdata", "exec", "repeatAlignment_bsmap.sh"), package="epiRepeatR")
						curExec <- "bsmap"
						curOtherArgs <- .config$alignment.params.bs
					} else if (alnType == "chip_bwa"){
						curScript <- system.file(file.path("extdata", "exec", "repeatAlignment_bwa_chip.sh"), package="epiRepeatR")
						curExec <- "bwa"
						curOtherArgs <- .config$alignment.params.chip
					} else {
						logger.warning(c("Did not find aligner for step", stepId,"--> skipping alignment step"))
					}
					doStep <- inputPresent & !is.null(curScript)
					if (doStep){
						if (alnType == "chip_bwa"){
							do.bwaIndex <- TRUE
						}
						if (sum(stepInds.input)>1){
							logger.warning(c("Multiple input bam files found for step",stepId,"--> picking the first one"))
						}
						inBamFn  <- ft[which(stepInds.input)[1],"fileName"]
						outBamFn <- file.path(paste0("${STEPDIR:", stepName, "}"), paste0(stepId,".bam"))


						curArgs <- unname(c(
							curScript,
							curExec, .config$samtools.exec,
							inBamFn, outBamFn,
							.config$refFasta,
							file.path("${TEMPDIR}", stepId),
							curOtherArgs
						))

						args.repeatAlignment <- c(args.repeatAlignment, list(curArgs))

						ft <- addFile(ft, outBamFn, dn, mn, sn, stepName)
					}
				}
			}
		}
		if (length(cmd.repeatAlignment) > 0 && length(args.repeatAlignment) > 0){
			# check if the reference has already been indexed for BWA alignment (ChIP) and if not, do it now
			if (do.bwaIndex) indexRepeatReference.bwa()

			pipr <- addStep(pipr, stepName, cmd.repeatAlignment, args.repeatAlignment, parents="bamExtract")
		}
		#-----------------------------------------------------------------------
		stepName <- "methCalling"
		logger.status(c("step:", stepName))
		cmd.methCalling <- .config$rscript.exec
		args.methCalling <- list()
		rscript <- system.file(file.path("extdata", "exec", "methCalling.R"), package="epiRepeatR")
		for (sn in sampleNames){
			for (mn in sampleMarks[[sn]]){
				for (dn in getSampleMarkDatatypes(.Object, sn, mn)){
					stepId <- paste(sn,mn,dn,stepName,sep="_")
					expInds <- ft[,"sampleName"]==sn & ft[,"mark"]==mn & ft[,"dataType"]==dn
					stepInds.input <- expInds & ft[,"fileType"]=="bam" & (ft[,"analysisStep"] %in% c("repeatAlignment")) & is.element(dn,c("WGBS","RRBS"))
					inputPresent <- any(stepInds.input)
					doStep <- inputPresent
					if (doStep){
						if (sum(stepInds.input)>1){
							logger.warning(c("Multiple input bam files found for step",stepId,"--> picking the first one"))
						}
						inBamFn  <- ft[which(stepInds.input)[1],"fileName"]
						outFn <- file.path(paste0("${STEPDIR:", stepName, "}"), paste0(stepId,".rds"))

						curArgs <- unname(c(
							rscript,
							"--in", inBamFn,
							"--out", outFn,
							"--config", cfgSavePath.forPipe
						))
						args.methCalling <- c(args.methCalling, list(curArgs))

						ft <- addFile(ft, outFn, dn, mn, sn, stepName)
					}
				}
			}
		}
		if (length(cmd.methCalling) > 0 && length(args.methCalling) > 0){		
			pipr <- addStep(pipr, stepName, cmd.methCalling, args.methCalling, parents="repeatAlignment")
		}
		#-----------------------------------------------------------------------
		stepName <- "chipQuantification"
		logger.status(c("step:", stepName))
		cmd.chipQuantification <- .config$rscript.exec
		args.chipQuantification <- list()
		rscript <- system.file(file.path("extdata", "exec", "chipQuantification.R"), package="epiRepeatR")
		for (sn in sampleNames){
			marks4sample <- sort(unique(ft[ft[, "sampleName"]==sn & ft[, "dataType"]=="ChIPseq", "mark"]))
			inds.input <- ft[,"sampleName"]==sn & ft[, "dataType"]=="Input" & (ft[, "analysisStep"] %in% c("repeatAlignment"))
			bamInput <- ft[which(inds.input)[1], "fileName"]
			chipBams <- sapply(marks4sample ,FUN=function(mm){
				idx <- ft[,"sampleName"]==sn & ft[,"dataType"]=="ChIPseq" & ft[,"mark"]==mm & (ft[,"analysisStep"] %in% c("repeatAlignment"))
				return(ft[which(idx)[1],"fileName"])
			})
			names(chipBams) <- marks4sample
			doStep <- any(inds.input) && length(chipBams) > 0
			if (doStep){
				for (mn in marks4sample){
					stepId <- paste(sn, mn, stepName, sep="_")
					outFn <- file.path(paste0("${STEPDIR:", stepName, "}"), paste0(stepId, ".rds"))
					curArgs <- unname(c(
						rscript,
						"--in", bamInput,
						"--chip", chipBams[mn],
						"--out", outFn,
						"--config", cfgSavePath.forPipe
					))
					args.chipQuantification <- c(args.chipQuantification, list(curArgs))
					ft <- addFile(ft, outFn, "ChIPseq", mn, sn, stepName)
				}
			}
		}
		if (length(cmd.chipQuantification) > 0 && length(args.chipQuantification) > 0){
			pipr <- addStep(pipr, stepName, cmd.chipQuantification, args.chipQuantification, parents="repeatAlignment")
		}
		#-----------------------------------------------------------------------
		stepName <- "repeatAlignmentStats"
		logger.status(c("step:", stepName))
		cmd.repeatAlignmentStats <- .config$rscript.exec
		rscript <- system.file(file.path("extdata", "exec", "repeatAlignmentStats.R"), package="epiRepeatR")
		inputPresent <- (ft[,"analysisStep"] %in% c("bamExtract", "repeatAlignment")) & ft[,"fileType"]=="bam"
		doStep <- any(inputPresent)
		if (doStep){
			outDir <- file.path(paste0("${STEPDIR:", stepName, "}"))
			ft.sub <- ft[inputPresent,]
			# ft.sub[,"SMT_ID"] <- paste(ft.sub[,"sampleName"], ft.sub[,"mark"], ft.sub[,"dataType"], sep="_") #an id based on sampleName, mark, dataType
			ft.be  <- ft.sub[ft.sub[,"analysisStep"] == "bamExtract",]
			ft.ra  <- ft.sub[ft.sub[,"analysisStep"] == "repeatAlignment",]

			ft.join <- merge(ft.be, ft.ra, by=c("dataType", "mark", "sampleName"), suffixes = c(".bamExtract", ".repeatAlignment"))

			inFileTable <- ft.join[, c("dataType", "mark", "sampleName", "fileName.bamExtract", "fileName.repeatAlignment")]
			inFileTable[,"fileName.bamExtract"]      <- muPipeR:::parseJobStrings(pipr, inFileTable[,"fileName.bamExtract"])
			inFileTable[,"fileName.repeatAlignment"] <- muPipeR:::parseJobStrings(pipr, inFileTable[,"fileName.repeatAlignment"])
			# TODO: avoid full file names in input file table
			inFileTable.fn <- paste0(stepName, "_inputFiles.tsv")
			inFileTable.fn.forPipe <- file.path(cfgDir.forPipe, inFileTable.fn)
			inFileTable.fn         <- file.path(cfgDir, inFileTable.fn)
			write.table(inFileTable, file=inFileTable.fn, quote=FALSE, row.names=FALSE, sep="\t", col.names=TRUE)

			args.repeatAlignmentStats <- c(
				rscript,
				"--in", inFileTable.fn.forPipe,
				"--out", outDir,
				"--config", cfgSavePath.forPipe
			)

			pipr <- addStep(pipr, stepName, cmd.repeatAlignmentStats, args.repeatAlignmentStats, parents=c("bamExtract", "repeatAlignment"))
		}
		#-----------------------------------------------------------------------
		stepName <- "plotRepeatGroupTreesMeth"
		logger.status(c("step:", stepName))
		cmd.plotRepeatGroupTreesMeth <- .config$rscript.exec
		rscript <- system.file(file.path("extdata", "exec", "plotRepeatGroupTreesMeth.R"), package="epiRepeatR")
		inputPresent <- (ft[,"mark"] %in% c("DNAmeth")) & (ft[,"analysisStep"] %in% c("methCalling")) & ft[,"fileType"]=="rds"
		doStep <- any(inputPresent)
		if (doStep){
			outDir <- file.path(paste0("${STEPDIR:", stepName, "}"))
			sampleNames <- ft[inputPresent, "sampleName"]
			if (any(duplicated(sampleNames))){
				logger.error("Multiple methylation call files found for a number of samples")
			}
			inFns <- ft[inputPresent,"fileName"]
			inFns.full <- muPipeR:::parseJobStrings(pipr, inFns)
			# TODO: avoid full file names in input file table
			inFileTable <- data.frame(fileName=inFns.full, sampleName=sampleNames, stringsAsFactors=FALSE)
			inFileTable.fn <- paste0(stepName, "_inputFiles.tsv")
			inFileTable.fn.forPipe <- file.path(cfgDir.forPipe, inFileTable.fn)
			inFileTable.fn         <- file.path(cfgDir, inFileTable.fn)
			write.table(inFileTable, file=inFileTable.fn, quote=FALSE, row.names=FALSE, sep="\t", col.names=TRUE)

			args.plotRepeatGroupTreesMeth <- c(
				rscript,
				"--in", inFileTable.fn.forPipe,
				"--out", outDir,
				"--config", cfgSavePath.forPipe,
				"--anaman", anamanSavePath.forPipe
			)

			pipr <- addStep(pipr, stepName, cmd.plotRepeatGroupTreesMeth, args.plotRepeatGroupTreesMeth, parents="methCalling")
		}
		#-----------------------------------------------------------------------
		stepName <- "plotRepeatMarkTree"
		logger.status(c("step:", stepName))
		cmd.plotRepeatMarkTree <- .config$rscript.exec
		rscript <- system.file(file.path("extdata", "exec", "plotRepeatMarkTree.R"), package="epiRepeatR")
		inputPresent <- (ft[,"analysisStep"] %in% c("methCalling", "chipQuantification")) & ft[,"fileType"]=="rds"
		# carry out the step if there are chipQuantifications
		doStep <- any(inputPresent) && any(ft[inputPresent,"analysisStep"]=="chipQuantification")
		if (doStep){
			outDir <- file.path(paste0("${STEPDIR:", stepName, "}"))
			sampleNames <- ft[inputPresent, "sampleName"]
			markNames   <- ft[inputPresent,"mark"]
			sampleNames.ext <- paste(sampleNames, markNames, sep="_")
			if (any(duplicated(sampleNames.ext))){
				logger.error("Multiple quantification files found for the same sample_mark combination present")
			}
			quantFns <- ft[inputPresent,"fileName"]
			inFns.full <- muPipeR:::parseJobStrings(pipr, quantFns)
			# TODO: avoid full file names in input file table
			inFileTable <- data.frame(fileName=inFns.full, sampleName=sampleNames, markName=markNames, stringsAsFactors=FALSE)
			inFileTable.fn <- paste0(stepName, "_inputFiles.tsv")
			inFileTable.fn.forPipe <- file.path(cfgDir.forPipe, inFileTable.fn)
			inFileTable.fn         <- file.path(cfgDir, inFileTable.fn)
			write.table(inFileTable, file=inFileTable.fn, quote=FALSE, row.names=FALSE, sep="\t", col.names=TRUE)

			args.plotRepeatMarkTree <- c(
				rscript,
				"--in", inFileTable.fn.forPipe,
				"--out", outDir,
				"--config", cfgSavePath.forPipe,
				"--anaman", anamanSavePath.forPipe
			)

			pipr <- addStep(pipr, stepName, cmd.plotRepeatMarkTree, args.plotRepeatMarkTree, parents=c("methCalling", "chipQuantification"))
		}

		write.table(ft, file=fileTabPath, quote=FALSE, row.names=FALSE, sep="\t", col.names=TRUE)
		logger.info(c("Saved file table to", fileTabPath))
		
		return(pipr)
	}
)

#' readFileTable
#'
#' Read a file table and return the corresponding AnalysisManager class object
#'
#' @param fileName	Filename of the Analysis configuration table. It should contain columns for filenames of the input files (column name: 'fileName'), the data type of the input files (column name: 'dataType'), a sample identifier ('sampleName'), an epigenetic mark identifier column ('mark') and a column indicating the analysis Step that generated the input file ('analysisStep')
#' @param groupColumns	column names of sample annotation. Can be set to NULL to use all columns not in the required input file annotation columns
#' @return a \code{\linkS4class{AnalysisManager}} object describing the analysis
#'
#' @author Fabian Mueller
#' @noRd
readFileTable <- function(fName, groupColumns=c()) {
	columnMap=c(fileName="fileName", dataType="dataType", sampleName="sampleName", mark="mark", analysisStep="analysisStep")

	tstTab <-read.table(fName, sep="\t", header=TRUE, stringsAsFactors=FALSE,nrows=1)
	if (is.null(groupColumns)){
		groupColumns <- setdiff(colnames(tstTab), columnMap)
	}

	new("AnalysisManager", fName, groupColumns, columnMap)
}

#' getSampleGroups
#'
#' given a sample annotation table, return the sample groups
#'
#' @param anno	sample annotation table. must contain sample ids as row names
#' @param addAll	flag indicating whether a group containing all samples in the annotation table should also be created
#' @param addIndividual	flag indicating whether a grouping containing a seperate group for each individual sample should be added
#' @return a list containing lists of sample identifiers. The first level list contains lists for each applicable column in the annotation table. On the second level, each list contains the sample ids for each group identified in that column 
#'
#' @details
#' gets the annotation from all columns which don't have only unique or only identical values
#'
#' @author Fabian Mueller
#' @noRd
getSampleGroups <- function(anno, addAll=FALSE, addIndividual=FALSE) {
	res <- list()
	for (j in 1:ncol(anno)){
		ccn <- colnames(anno)[j]
		vv <- anno[,j]
		isValidGrpCol <- (length(unique(vv)) < length(vv)) && (length(unique(vv)) > 1)
		if (isValidGrpCol){
			grps <- lapply(unique(vv), FUN=function(lvl){
				rownames(anno)[vv==lvl]
			})
			names(grps) <- unique(vv)

			res_cur <- list()
			res_cur[[ccn]] <- grps
			res <- c(res, res_cur)
		}
	}
	if (addAll){
		res_cur <- list()
		res_cur[["ALL"]] <- list("ALL"=rownames(anno))
		res <- c(res, res_cur)
	}
	if (addIndividual){
		res_cur <- lapply(rownames(anno), identity)
		names(res_cur) <- rownames(anno)
		res_cur <- list(".sample"=res_cur)
		res <- c(res, res_cur)
	}
	return(res)
}
