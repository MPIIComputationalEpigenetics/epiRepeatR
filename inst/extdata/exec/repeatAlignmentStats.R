library(argparse)
suppressPackageStartupMessages(library(epiRepeatR))
ap <- ArgumentParser()
ap$add_argument("-i", "--in", action="store", dest="inFileTable", help="Input file table (tab-separated).")
ap$add_argument("-o", "--out", action="store", dest="output", help="Output directory")
ap$add_argument("-c", "--config", action="store", help="Config file (json)")
cmdArgs <- ap$parse_args() #problem: too long of a command line

# saveRDS(cmdArgs, file.path(cmdArgs$output, "cmdargs.rds"))

outDir <- cmdArgs$output
loadConfig(cmdArgs$config)
inFileTable <- read.table(cmdArgs$inFileTable, sep="\t", comment.char="", header=TRUE, stringsAsFactors=FALSE)

logger.start("Getting read counts")
	alnStats <- do.call("rbind", lapply(1:nrow(inFileTable), function(i){
		epiRepeatR:::getAlnStats(inFileTable[i, "fileName.repeatAlignment"], inFileTable[i, "fileName.bamExtract"])
	}))
	alnStats <- data.frame(inFileTable[,c("sampleName", "mark", "dataType")], alnFile=inFileTable[,"fileName.repeatAlignment"],  alnStats, stringsAsFactors=FALSE)
logger.completed()

logger.start("Writing output")
	fn <- file.path(outDir, "alignmentStats.tsv")
	write.table(alnStats, file=fn, quote=FALSE, row.names=FALSE, sep="\t", col.names=TRUE)
logger.completed()
