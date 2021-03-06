library(argparse)
suppressPackageStartupMessages(library(epiRepeatR))
library(tools)
ap <- ArgumentParser()
ap$add_argument("-i", "--in", action="store", dest="inFileTable", help="Input file table (tab-separated).")
ap$add_argument("-o", "--out", action="store", dest="output", help="Output file name")
ap$add_argument("-c", "--config", action="store", help="Config file (json)")
cmdArgs <- ap$parse_args() #problem: too long of a command line
logger.cmd.args(cmdArgs)

loadConfig(cmdArgs$config)

if (is.element("debug", names(epiRepeatR:::.config)) && epiRepeatR:::.config$debug){
	saveRDS(cmdArgs, file.path(dirname(cmdArgs$output), "cmdargs.rds"))
}

outBam <- cmdArgs$output

fns <- read.table(cmdArgs$inFileTable, sep="\t", comment.char="", header=FALSE, stringsAsFactors=FALSE)[,1]

if (all(file_ext(fns)=="bam")){
	samtoolsExec <- epiRepeatR:::.config$samtools.exec
	args <- c("merge", "-r", outBam, fns)

	logger.start("Merging bam files")
		res <- system2(samtoolsExec, args)
	logger.completed()

	logger.start("Indexing output file")
		res.index <- system2(samtoolsExec, c("index", outBam))
	logger.completed()
} else if (all(grepl("\\.(fq|fastq)(\\.gz)?$", fns))){
	logger.start("Merging fastq files")
		args <- c(fns, ">", outBam)
		res <- system2("cat", args)
	logger.completed()
} else {
	logger.error("Mixed input file types are not supported")
}
