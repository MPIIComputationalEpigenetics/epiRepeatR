library(argparse)
suppressPackageStartupMessages(library(epiRepeatR))
ap <- ArgumentParser()
ap$add_argument("-i", "--in", action="store", dest="input", help="Input file (RDS) containing a RepeatEpigenomeCollection object as R dataset.")
ap$add_argument("-o", "--out", action="store", dest="output", help="Output directory")
ap$add_argument("-c", "--config", action="store", help="Config file (json)")
cmdArgs <- ap$parse_args() #problem: too long of a command line
logger.cmd.args(cmdArgs)

loadConfig(cmdArgs$config)

if (is.element("debug", names(epiRepeatR:::.config)) && epiRepeatR:::.config$debug){
	saveRDS(cmdArgs, file.path(cmdArgs$output, "cmdargs.rds"))
}

rec <- readRDS(cmdArgs$input)

epiRepeatR:::createRepPlot_markTree(
	rec, cmdArgs$output,
	dendroMethod=epiRepeatR:::.config$plotRepTree.dendroMethod,
	minReads=epiRepeatR:::.config$plotRepTree.meth.minReads, minCpGs=epiRepeatR:::.config$plotRepTree.meth.minCpGs
)
