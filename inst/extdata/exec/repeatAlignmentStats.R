library(argparse)
suppressPackageStartupMessages(library(epiRepeatR))
library(ggplot2)
ap <- ArgumentParser()
ap$add_argument("-i", "--in", action="store", dest="inFileTable", help="Input file table (tab-separated).")
ap$add_argument("-o", "--out", action="store", dest="output", help="Output directory")
ap$add_argument("-c", "--config", action="store", help="Config file (json)")
cmdArgs <- ap$parse_args() #problem: too long of a command line
logger.cmd.args(cmdArgs)

loadConfig(cmdArgs$config)

if (is.element("debug", names(epiRepeatR:::.config)) && epiRepeatR:::.config$debug){
	saveRDS(cmdArgs, file.path(cmdArgs$output, "cmdargs.rds"))
}

outDir <- cmdArgs$output
inFileTable <- read.table(cmdArgs$inFileTable, sep="\t", comment.char="", header=TRUE, stringsAsFactors=FALSE)

logger.start("Getting read counts and stats")
	alnStats <- do.call("rbind", lapply(1:nrow(inFileTable), function(i){
		logger.status(c("Processing (",i,"):", paste(inFileTable[i, c("sampleName","mark","dataType")], collapse=" - ")))
		return(epiRepeatR:::getAlnStats(inFileTable[i, "fileName.repeatAlignment"], inFileTable[i, "fileName.bamExtract"]))
	}))
	alnStats <- data.frame(inFileTable[,c("sampleName", "mark", "dataType")], alnFile=inFileTable[,"fileName.repeatAlignment"],  alnStats, stringsAsFactors=FALSE)
logger.completed()

logger.start("Writing output")
	fn <- file.path(outDir, "alignmentStats.tsv")
	write.table(alnStats, file=fn, quote=FALSE, row.names=FALSE, sep="\t", col.names=TRUE)
logger.completed()

logger.start("Generating Plots")
	theme_set(theme_bw())
	fn <- file.path(outDir, "alignmentStats.pdf")
	pp <- plotRepeatAlignmentStats(alnStats)
	ggsave(fn, pp, width=10,height=10)
logger.completed()
