library(argparse)
suppressPackageStartupMessages(library(epiRepeatR))
ap <- ArgumentParser()
ap$add_argument("-i", "--in", action="store", dest="input", help="Input file (bam)")
ap$add_argument("-o", "--out", action="store", dest="output", help="Output file (rds)")
ap$add_argument("-c", "--config", action="store", help="Config file (json)")
cmdArgs <- ap$parse_args()

loadConfig(cmdArgs$config)

ra <- epiRepeatR:::RepeatAlignmentBiSeq(cmdArgs$input)
methCallObj <- epiRepeatR:::getMethylationCalls(ra)
saveRDS(methCallObj, cmdArgs$output)
