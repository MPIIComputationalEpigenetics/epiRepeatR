library(argparse)
suppressPackageStartupMessages(library(epiRepeatR))
ap <- ArgumentParser()
ap$add_argument("-i", "--in", action="store", dest="inFileTable", help="Input file table (tab-separated). Should contain three columns (with headers): fileName, sampleName, markName")
ap$add_argument("-o", "--out", action="store", dest="output", help="Output directory")
ap$add_argument("-c", "--config", action="store", help="Config file (json)")
ap$add_argument("-a", "--anaman", action="store", help="Analysis manager object (rds)")
cmdArgs <- ap$parse_args() #problem: too long of a command line

# saveRDS(cmdArgs, file.path(cmdArgs$output, "cmdargs.rds"))

loadConfig(cmdArgs$config)
anaMan <- readRDS(cmdArgs$anaman)
inFileTable <- read.table(cmdArgs$inFileTable, sep="\t", comment.char="", header=TRUE, stringsAsFactors=FALSE)
inFiles <- inFileTable[,"fileName"]
sampleNames <- inFileTable[,"sampleName"]
markNames <- inFileTable[,"markName"]

sampleTab <- epiRepeatR:::getSampleAnnot(anaMan)
invalidSamples <- setdiff(sampleNames,rownames(sampleTab))
if (length(invalidSamples)>0){
	logger.error(c("The following samples do not exist in the annotation table:",paste(invalidSamples,collapse=",")))
}
sampleTab <- sampleTab[sampleNames,,drop=FALSE]
ggs <- epiRepeatR:::getSampleGroups(sampleTab, addAll=TRUE)

epiRepeatR:::createRepPlot_markTree(
	inFiles, sampleNames, markNames, cmdArgs$output,
	sampleGroups=ggs,
	dendroMethod=epiRepeatR:::.config$plotRepTree.dendroMethod,
	minReads=epiRepeatR:::.config$plotRepTree.meth.minReads, minCpGs=epiRepeatR:::.config$plotRepTree.meth.minCpGs
)
