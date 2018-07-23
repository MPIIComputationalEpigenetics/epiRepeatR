.config <- new.env()
.config$n.processes <- 0
.config$rscript.exec <- "Rscript"
.config$species <-  "human"
.config$assembly <-  NULL
# .config$refFasta <- "/TL/deep-external01/archive00/references/assemblies/RepBase19.07/human.fasta/humrep.fa"
.config$refFasta <- "/TL/deep-external01/archive00/references/assemblies/RepBase20.04/Homo_sapiens_all.fa"
.config$refEmbl  <- NULL 
.config$aligner.bs <- "bsmap"
.config$alignment.params.bs <- "-g 3 -v 0.2"
.config$aligner.chip <- "chip_bwa" # ChIP-seq aligner
.config$alignment.params.chip.bwa <- "-t 8 -q 20 -b" # to be appended to BWA ChIP-seq aligner call
.config$alignment.params.chip.bowtie2 <- "--local --sensitive-local -p 8" # to be appended to Bowtie2 ChIP-seq aligner call
.config$aligner.atac <- "atac_bt2" # ATAC-seq aligner
.config$alignment.params.atac.bowtie2 <- "--local --sensitive-local -p 8" # to be appended to Bowtie2 ChIP-seq aligner call
.config$samtools.exec <- "samtools" # location of the samtools executable
.config$chip.mergeInput <- FALSE
.config$tempDir <- tempdir()
.config$inputBam.mappingStatus <- "all"
.config$plotRepTree.dendroMethod <- "repeatFamily"
.config$plotRepTree.normEnrich <- "none"
.config$plotRepTree.meth.minCpGs  <- 2
.config$plotRepTree.meth.minReads <- 100
.config$annotCols.replicates <- NULL
.config$meth.minCpGcov  <- 5 #NULL
.config$meth.methCallFormat  <- "BisSNP" #"EPP"
.config$genomeRepeatTrack <- NULL
.config$debug <- FALSE
#assignInNamespace(".config", .config, "epiRepeatR") # to manually assign in non-exported variables

#' setConfigElement
#'
#' Set a configuration item to a given value
#'
#' @param name	name of the config item
#' @param value	value of the config item
#' @return nothing of particular interest.
#' @section Options used by the package:
#' \describe{
#'   \item{\bold{\code{n.processes}}\code{ = 0}}{
#'        Number of processes used during parallel computations.
#'   }
#'   \item{\bold{\code{rscript.exec}}\code{ = "Rscript"}}{
#'        \code{Rscript executable.}
#'   }
#'   \item{\bold{\code{species}}\code{ = "human"}}{
#'        The organism/taxon the analysis pertains to.
#'   }
#'   \item{\bold{\code{refFasta}}\code{ = } MPII default location}{
#'        The reference fasta file for repeat elements. See vignette for instructions on how they can be obtained.
#'   }
#'   \item{\bold{\code{aligner.bs}}\code{ = "bsmap"}}{
#'        Aligner used for mapping bisulfite sequencing reads to the repeat reference. Currently only "bsmap" is supported.
#'   }
#'   \item{\bold{\code{alignment.params.bs}}\code{ = "-g 3 -v 0.2"}}{
#'        Additional parameters appended during the command line call to the bisulfite aligner.
#'   }
#'   \item{\bold{\code{aligner.chip}}\code{ = "bwa_chip"}}{
#'        Aligner used for mapping ChIP-seq reads to the repeat reference. Currently only "chip_bwa" is supported.
#'   }
#'   \item{\bold{\code{alignment.params.chip.bwa}}\code{ = "-t 8 -q 20 -b"}}{
#'        Additional parameters appended during the command line call to the BWA ChIP-seq aligner.
#'   }
#'   \item{\bold{\code{samtools.exec}}\code{ = "samtools"}}{
#'        Location of the samtools executable.
#'   }
#'   \item{\bold{\code{chip.mergeInput}}\code{ = FALSE}}{
#'        Flag indicating whether all bam files containing Input/WCE for ChIP-seq should be merged into one joint Input.
#'   }
#'   \item{\bold{\code{tempDir}}\code{ = tempdir()}}{
#'         Temporary directory for analysis. Specify to be the empty string ("") to create a temp directory in the analysis directory during runAnalysis.
#'   }
#'   \item{\bold{\code{inputBam.mappingStatus}}\code{ = "all"}}{
#'        Which reads should be extracted from the bam files and subsequently mapped to repetitive elements:
#'        Only reads also mapped in the source bam file ("mapped"), only reads unmapped in the source bam file ("unmapped"),
#'        or both ("all", default).
#'   }
#'   \item{\bold{\code{plotRepTree.dendroMethod}}\code{ = "repeatFamily"}}{
#'         Method for plotting the repeat subfamily dendrogram. Valid methods include 
#'                  grouping by repeat subfamily ("repeatFamily"),
#'                  grouping by hierarchical clustering based on k-mer counts in the repeat sequence (Euclidean distance, complete linkage) ("hierClust"),
#'                  grouping by hierarchical clustering based on occurrences of terms in the annotation fields of a repeat 
#'                  ("annotClust"; requires that the repeat references has been annotated from the EMBL format.)
#'   }
#'   \item{\bold{\code{plotRepTree.normEnrich}}\code{ = "none"}}{
#'        Method for normalizing enrichment data before plotting. Currently supported are:
#'               \code{none} (no normalization),
#'               \code{standard} (subtract the mean, devide by standard deviation),
#'               \code{scale} (scale to the interval [0,1]) and
#'               \code{quantile} (Quantile normalization)
#'   }
#'   \item{\bold{\code{plotRepTree.meth.minCpGs}}\code{ = 2}}{
#'        Minimum number of CpGs to be contained in a repeat element in order to be shown in the resulting methylation tree plots.
#'   }
#'   \item{\bold{\code{plotRepTree.meth.minReads}}\code{ = 100}}{
#'        Minimum number of reads that must match to a given repeat element in order to be shown in the resulting methylation tree plots.
#'   }
#'   \item{\bold{\code{annotCols.replicates}}\code{ = NULL}}{
#'        Column names or indices in the annotation column used for replicate analysis (in exploratory report)
#'   }
#'   \item{\bold{\code{meth.minCpGcov}}\code{ = 5}}{
#'        Minimum number of reads covering a CpG in a repeat element in order to be considered in computing methylation.
#'   }
#'   \item{\bold{\code{meth.methCallFormat}}\code{ = "BisSNP"}}{
#'        Format of input methylation calling files. Can be one of "BisSNP", "EPP"
#'   }
#'   \item{\bold{\code{genomeRepeatTrack}}\code{ = NULL}}{
#'      Path to an RDS file containing a GenomeRepeatTrack object. These files are used to map the sequencing reads instead of the consensus reference.
#'      Only meaningful if processing input files of type \code{genomeMethCalling} or \code{genomeAlignment}.
#'   }
#'   \item{\bold{\code{debug}}\code{ = FALSE}}{
#'        \code{Logical specifying whether the debug mode is enabled and additional debug-related output should be provided.}
#'   }
#' }
#' @author Fabian Mueller
#' @export
setConfigElement <- function(name, value){
	.config[[name]] <- value
}

#' getConfigElement
#'
#' Get the value for a configuration item
#'
#' @param name	name of the config item
#' @return the value of the config item
#' @author Fabian Mueller
#' @export
getConfigElement <- function(name){
	if (!exists(name, .config)){
		logger.warning(c("No such configuration element:", name, "--> NULL returned"))
	}
	.config[[name]]
}
#' saveConfig
#'
#' Save the current configuration to a configuration file (JSON)
#'
#' @param dest	Filename for the config file in JSON format.
#' @return nothing of particular interest.
#'
#' @author Fabian Mueller
#' @export
saveConfig <- function(dest){
	cat(toJSON(as.list(.config), pretty=TRUE), file=dest)
}

#' loadConfig
#'
#' Sets the configuration from a configuration file (JSON)
#'
#' @param cfgFile	Config file in JSON format. As output by \link{saveConfig}
#' @return nothing of particular interest. The configuration is set for the current environment
#'
#' @author Fabian Mueller
#' @export
loadConfig <- function(cfgFile){
	cfgList <- fromJSON(cfgFile)
	for (nn in names(cfgList)){
		if (is.element(nn,ls(.config))){
			.config[[nn]] <- cfgList[[nn]]
		} else {
			logger.warning(c("Ignoring unknown item '",nn,"' in when loading configuration from file",cfgFile))
		}
	}
}

#' setUpParallelFromConfig
#'
#' Inspect the current configuration and set up a parallel environment using \code{RnBeads::parallel.setup}
#'
#' @return nothing of particular interest.
#'
#' @author Fabian Mueller
#' @noRd
setUpParallelFromConfig <- function(){
	if (!logger.isinitialized()) logger.start(fname=NA)
	#set up multicore
	if (.config$n.processes > 0){
		parallel.setup(.config$n.processes)
	} else {
		logger.info(c("Using single core"))
	}
}
#' getArgParser
#'
#' Gets an argument parser for command line arguments for the repeat analysis pipeline
#'
#' @return an \code{ArgumentParser} (see \code{argparse} package) object
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' ap <- getArgParser()
#' ap$print_help()
#' }
getArgParser <- function(){
	require(argparse)
	ap <- ArgumentParser()
	# !!! the dest have to be consistent with the names of .config or the arguments of runAnalysis
	ap$add_argument("-a","--analysisDir", action="store", dest="anaDir", help="Directory where analysis will be performed and results will be stored. NOTE: the directory structure will be managed by epiRepeatR.")
	ap$add_argument("-f","--fileTable", action="store", dest="fileTable", help="File table containing paths and anotation for the input files in the analysis.")
	ap$add_argument("-c","--configFile", action="store", dest="cfgFile", help="JSON config file specifying further analysis options.")
	ap$add_argument("--resetToStep", action="store", dest="resetToStep", help="Reset the analysis to the beginning of this step prior to running the analysis.")
	ap$add_argument("-r","--referenceFasta", action="store", dest="refFasta", help="Fasta file containing the reference sequences for repeat elements.")
	ap$add_argument("-p","--numProcesses", action="store", dest="n.processes", type="integer", help="Number pf processes to be used for analysis.")	
	ap$add_argument("--alignerBS", action="store", dest="aligner.bs", help="Aligner to be used for mapping bisulfite reads to reference repeats. Currently only 'bsmap' is supported")
	ap$add_argument("--alignerChip", action="store", dest="aligner.chip", help="Aligner to be used for mapping ChIP-seq reads to reference repeats. Currently only 'chip_bwa' is supported")
	ap$add_argument("--alignerAtac", action="store", dest="aligner.atac", help="Aligner to be used for mapping ATAC-seq reads to reference repeats. Currently only 'atac_bwa' is supported")
	ap$add_argument("--inputBamUnmapped", action="store_true", dest="inputBamUnmapped", help="Only process unmapped reads from an input bam file")
	return(ap)
}
#' loadConfigFromCommandLine
#'
#' Sets the configuration from command line options
#'
#' @return a list of parameters for input to the \link{runAnalysis} function
#'
#' @author Fabian Mueller
#' @export
loadConfigFromCommandLine <- function(){
	ap <- getArgParser()
	cmdArgs <- ap$parse_args()
	# first look for a config file in JSON format
	if (!is.null(cmdArgs$cfgFile)){
		loadConfig(cfgFile)
	}
	# then set options from command line
	cfgItems <- intersect(names(cmdArgs),ls(.config))
	for (nn in cfgItems){
		curCfg <- cmdArgs[[nn]]
		if (!is.null(curCfg)){
			.config[[nn]] <- curCfg
		}		
	}
	#process options
	if (cmdArgs$inputBamUnmapped){
		.config$inputBam.mappingStatus <- "unmapped"
	}
	#delete options unrecognized by runAnalysis(...)
	cmdArgs[["inputBamUnmapped"]] <- NULL
	cmdArgs[["cfgFile"]] <- NULL

	otherItems <- setdiff(names(cmdArgs),ls(.config))
	return(cmdArgs[otherItems])
}
