library(argparse)
suppressPackageStartupMessages(library(epiRepeatR))
ap <- ArgumentParser()
ap$add_argument("-i", "--in", action="store", dest="inFileTable", help="Input file table (tab-separated). Should contain three columns (with headers): fileName, sampleName, markName")
ap$add_argument("-o", "--out", action="store", dest="output", help="Output directory")
ap$add_argument("-c", "--config", action="store", help="Config file (json)")
cmdArgs <- ap$parse_args() #problem: too long of a command line

# saveRDS(cmdArgs, file.path(cmdArgs$output, "cmdargs.rds"))

loadConfig(cmdArgs$config)
inFileTable <- read.table(cmdArgs$inFileTable, sep="\t", comment.char="", header=TRUE, stringsAsFactors=FALSE)
inFiles <- inFileTable[,"fileName"]
sampleNames <- inFileTable[,"sampleName"]
markNames <- inFileTable[,"markName"]


epiRepeatR:::createRepPlot_markTree(
	inFiles, sampleNames, markNames, cmdArgs$output,
	minReads=epiRepeatR:::.config$plotRepTree.meth.minReads, minCpGs=epiRepeatR:::.config$plotRepTree.meth.minCpGs
)
