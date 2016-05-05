library(argparse)
suppressPackageStartupMessages(library(epiRepeatR))
ap <- ArgumentParser()
ap$add_argument("-i", "--in", action="store", dest="input", help="Input file (bam;  bed for genome meth calling files)")
ap$add_argument("-o", "--out", action="store", dest="output", help="Output file (rds)")
ap$add_argument("-c", "--config", action="store", help="Config file (json)")
ap$add_argument("-g", "--grt", action="store", dest="genomeRepTrack", help="Enables parsing from genome methylation calls. Specifies Path to an RDS file containing a GenomeRepeatTrack object.")
cmdArgs <- ap$parse_args()
logger.cmd.args(cmdArgs)

loadConfig(cmdArgs$config)

if (!is.null(cmdArgs$genomeRepTrack)){
	logger.info("Calling methylation from genome methylation calls")
	grt <- readRDS(cmdArgs$genomeRepTrack)
	gmc <- epiRepeatR:::GenomeMethylationCalls(cmdArgs$input, grt)
	fn <- paste0(gsub("\\.rds$", "", cmdArgs$output, ignore.case=TRUE), "_gmc.rds")
	saveRDS(gmc, fn)
	methCallObj <- epiRepeatR:::getMethylationCalls(gmc)
} else {
	logger.info("Calling methylation from repeat alignment")
	ra <- epiRepeatR:::RepeatAlignmentBiSeq(cmdArgs$input)
	methCallObj <- epiRepeatR:::getMethylationCalls(ra)
}
saveRDS(methCallObj, cmdArgs$output)
