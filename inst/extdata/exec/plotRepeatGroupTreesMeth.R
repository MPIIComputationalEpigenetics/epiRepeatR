library(argparse)
suppressPackageStartupMessages(library(epiRepeatR))
ap <- ArgumentParser()
ap$add_argument("-i", "--in", action="store", dest="input", help="Input file (RDS) containing a RepeatEpigenomeCollection object as R dataset.")
ap$add_argument("-o", "--out", dest="output", action="store", help="Output Prefix")
ap$add_argument("-c", "--config", action="store", help="Config file (json)")
cmdArgs <- ap$parse_args()
logger.cmd.args(cmdArgs)

loadConfig(cmdArgs$config)

if (is.element("debug", names(epiRepeatR:::.config)) && epiRepeatR:::.config$debug){
	saveRDS(cmdArgs, file.path(cmdArgs$output, "cmdargs.rds"))
}

rec <- readRDS(cmdArgs$input)

logger.info(c("using dendrogram method:",epiRepeatR:::.config$plotRepTree.dendroMethod))
epiRepeatR:::createRepPlot_groupSummaryTrees_meth(
	rec, cmdArgs$output,
	dendroMethod=getConfigElement("plotRepTree.dendroMethod"),
	minReads=getConfigElement("plotRepTree.meth.minReads"), minCpGs=getConfigElement("plotRepTree.meth.minCpGs"), minCpGcov=getConfigElement("meth.minCpGcov")
)
