library(argparse)
suppressPackageStartupMessages(library(epiRepeatR))
ap <- ArgumentParser()
ap$add_argument("-i", "--in", action="store", dest="input", help="Input file (bam)")
ap$add_argument("-j", "--chip", action="store", help="ChIP file (bam)")
ap$add_argument("-o", "--out", action="store", dest="output", help="Output file (rds)")
ap$add_argument("-c", "--config", action="store", help="Config file (json)")
ap$add_argument("-g", "--genome", action="store_true", dest="doGenomeAln", default=FALSE, help="Enables parsing from genome alignments.")
cmdArgs <- ap$parse_args()
logger.cmd.args(cmdArgs)

loadConfig(cmdArgs$config)

if (cmdArgs$doGenomeAln){
	logger.info("Quantify enrichment from genome alignment")
	ga.input <- readRDS(cmdArgs$input)
	ga.chip  <- readRDS(cmdArgs$chip)
	quantObj <- epiRepeatR:::computeEnrichment(ga.chip, ga.input)
} else {
	logger.info("Quantify enrichment from repeat alignment")
	ra.input <- epiRepeatR:::RepeatAlignment(cmdArgs$input)
	ra.chip  <- epiRepeatR:::RepeatAlignmentChip(cmdArgs$chip)
	quantObj <- epiRepeatR:::computeEnrichment(ra.chip, ra.input, useIdxStats=TRUE)
}
saveRDS(quantObj, cmdArgs$output)


