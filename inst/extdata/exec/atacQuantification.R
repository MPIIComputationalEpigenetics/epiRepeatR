library(argparse)
suppressPackageStartupMessages(library(epiRepeatR))
ap <- ArgumentParser()
ap$add_argument("-i", "--in", action="store", dest="input", help="Input file of aligned ATAC-seq reads (bam)")
ap$add_argument("-o", "--out", action="store", dest="output", help="Output file (rds)")
ap$add_argument("-c", "--config", action="store", help="Config file (json)")
ap$add_argument("-g", "--genome", action="store_true", dest="genomeRepTrack", default=FALSE, help="Enables parsing from genome methylation calls.")
cmdArgs <- ap$parse_args()
logger.cmd.args(cmdArgs)

loadConfig(cmdArgs$config)

if (cmdArgs$genomeRepTrack){
	logger.info("Quantify enrichment from genome alignment")
	ga.atac <- readRDS(cmdArgs$input)
	quantObj <- epiRepeatR:::normCounts(ga.atac)
} else {
	logger.info("Quantify enrichment from repeat alignment")
	ra.atac <- epiRepeatR:::RepeatAlignment(cmdArgs$input)
	quantObj <- epiRepeatR:::normCounts(ra.atac, useIdxStats=TRUE)
}
saveRDS(quantObj, cmdArgs$output)


