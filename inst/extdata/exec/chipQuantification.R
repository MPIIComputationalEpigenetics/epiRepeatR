library(argparse)
suppressPackageStartupMessages(library(epiRepeatR))
ap <- ArgumentParser()
ap$add_argument("-i", "--in", action="store", dest="input", help="Input file (bam)")
ap$add_argument("-j", "--chip", action="store", help="ChIP file (bam)")
ap$add_argument("-o", "--out", action="store", dest="output", help="Output file (rds)")
ap$add_argument("-c", "--config", action="store", help="Config file (json)")
cmdArgs <- ap$parse_args()

loadConfig(cmdArgs$config)

ra.input <- epiRepeatR:::RepeatAlignment(cmdArgs$input)
ra.chip  <- epiRepeatR:::RepeatAlignmentChip(cmdArgs$chip)
quantObj <- epiRepeatR:::computeEnrichment(ra.chip, ra.input, useIdxStats=TRUE)
saveRDS(quantObj, cmdArgs$output)


