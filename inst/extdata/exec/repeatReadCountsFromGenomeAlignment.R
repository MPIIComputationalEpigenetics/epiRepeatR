library(argparse)
suppressPackageStartupMessages(library(epiRepeatR))
ap <- ArgumentParser()
ap$add_argument("-i", "--in", action="store", dest="input", help="Input file (bam)")
ap$add_argument("-o", "--out", action="store", dest="output", help="Output file (rds)")
ap$add_argument("-c", "--config", action="store", help="Config file (json)")
ap$add_argument("-g", "--grt", action="store", dest="genomeRepTrack", help="[Optional] Enables parsing from genome methylation calls. Specifies Path to an RDS file containing a GenomeRepeatTrack object.")
cmdArgs <- ap$parse_args()
logger.cmd.args(cmdArgs)

loadConfig(cmdArgs$config)

if (is.null(cmdArgs$genomeRepTrack) || !file.exists(cmdArgs$genomeRepTrack)){
	logger.info(c("Creating genome repeat track for genome assembly", .config$assembly))
	grt <- GenomeRepeatTrack(.config$assembly)
} else {
	logger.info(c("Reading genome repeat track from", cmdArgs$genomeRepTrack))
	grt <- readRDS(cmdArgs$genomeRepTrack)
}
ga <- epiRepeatR:::GenomeAlignment(cmdArgs$input, repeatTrack=grt)
ga <- epiRepeatR:::storeReadCounts(ga)

saveRDS(ga, cmdArgs$output)


