suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(epiRepeatR))
ap <- ArgumentParser()
ap$add_argument("-i", "--in", action="store", dest="input", help="Input file of aligned ATAC-seq reads (bam)")
ap$add_argument("-o", "--out", action="store", dest="output", help="Output file (rds)")
ap$add_argument("-c", "--config", action="store", help="Config file (json)")
ap$add_argument("-g", "--genome", action="store_true", dest="doGenomeAln", default=FALSE, help="Enables parsing from genome alignments.")
ap$add_argument("--grt", action="store", dest="genomeRepTrack", help="Specifies Path to an RDS file containing a GenomeRepeatTrack object.")
cmdArgs <- ap$parse_args()
logger.cmd.args(cmdArgs)

# METHOD <- "zscore"
METHOD <- "genomeScale"

loadConfig(cmdArgs$config)

abund <- NULL
if (METHOD=="genomeScale"){
	if (!is.null(cmdArgs$genomeRepTrack)){
		logger.info("Retrieving abundances from genome repeat track")
		grt <- readRDS(cmdArgs$genomeRepTrack)
		abund <- epiRepeatR:::getRepeatGenomeCovg(grt)
	}
}
if (cmdArgs$doGenomeAln){
	logger.info("Quantify enrichment from genome alignment")
	ga.atac <- readRDS(cmdArgs$input)
	quantObj <- epiRepeatR:::normCounts(ga.atac, method=METHOD, abund=abund)
} else {
	logger.info("Quantify enrichment from repeat alignment")
	ra.atac <- epiRepeatR:::RepeatAlignment(cmdArgs$input)
	quantObj <- epiRepeatR:::normCounts(ra.atac, method=METHOD, abund=abund, useIdxStats=TRUE)
}
saveRDS(quantObj, cmdArgs$output)


