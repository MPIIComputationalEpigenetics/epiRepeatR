library(argparse)
suppressPackageStartupMessages(library(epiRepeatR))
ap <- ArgumentParser()
ap$add_argument("-i", "--in", action="store", dest="inFileTable", help="Input file table (tab-separated). Should contain three columns (with headers): fileName, sampleName, markName")
ap$add_argument("-o", "--out", action="store", dest="output", help="Output file name (RDS)")
ap$add_argument("-c", "--config", action="store", help="Config file (json)")
ap$add_argument("-a", "--anaman", action="store", help="Analysis manager object (rds)")
cmdArgs <- ap$parse_args() #problem: too long of a command line
logger.cmd.args(cmdArgs)

loadConfig(cmdArgs$config)

if (is.element("debug", names(epiRepeatR:::.config)) && epiRepeatR:::.config$debug){
	saveRDS(cmdArgs, file.path(dirname(cmdArgs$output), paste0(gsub(".rds$" ,"", basename(cmdArgs$output)), "_cmdargs.rds")))
}

anaMan <- readRDS(cmdArgs$anaman)
inFileTable <- read.table(cmdArgs$inFileTable, sep="\t", comment.char="", header=TRUE, stringsAsFactors=FALSE)
quantFns <- inFileTable[,"fileName"]
sampleNames <- inFileTable[,"sampleName"]
markNames <- inFileTable[,"markName"]

annot <- epiRepeatR:::getSampleAnnot(anaMan)
invalidSamples <- setdiff(sampleNames,rownames(annot))
if (length(invalidSamples)>0){
	logger.error(c("The following samples do not exist in the annotation table:",paste(invalidSamples,collapse=",")))
}

rec <- epiRepeatR:::RepeatEpigenomeCollection(quantFns, sampleNames, markNames, annot)
saveRDS(rec, cmdArgs$output)
