suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(epiRepeatR))
ap <- ArgumentParser()
ap$add_argument("-i", "--in", action="store", dest="input", help="Input file of aligned ATAC-seq reads (bam)")
ap$add_argument("-o", "--out", action="store", dest="output", help="Output file (rds)")
ap$add_argument("-c", "--config", action="store", help="Config file (json)")
ap$add_argument("-g", "--genome", action="store_true", dest="genomeRepTrack", default=FALSE, help="Enables parsing from genome methylation calls.")
ap$add_argument("-g", "--grt", action="store", dest="genomeRepTrack", help="Enables parsing from genome methylation calls. Specifies Path to an RDS file containing a GenomeRepeatTrack object.")
cmdArgs <- ap$parse_args()
logger.cmd.args(cmdArgs)

METHOD="zscore"

loadConfig(cmdArgs$config)

eCounts <- NULL
if (METHOD=="genomeScale"){
	if (!is.null(cmdArgs$genomeRepTrack)){
		grt <- readRDS(cmdArgs$genomeRepTrack)
		eCounts <- epiRepeatR:::getRepeatGenomeCovg(grt)
	}
}
if (cmdArgs$genomeRepTrack){
	logger.info("Quantify enrichment from genome alignment")
	ga.atac <- readRDS(cmdArgs$input)
	quantObj <- epiRepeatR:::normCounts(ga.atac, method=METHOD, eCounts=eCounts)
} else {
	logger.info("Quantify enrichment from repeat alignment")
	ra.atac <- epiRepeatR:::RepeatAlignment(cmdArgs$input)
	quantObj <- epiRepeatR:::normCounts(ra.atac, method=METHOD, eCounts=eCounts, useIdxStats=TRUE)
}
saveRDS(quantObj, cmdArgs$output)


