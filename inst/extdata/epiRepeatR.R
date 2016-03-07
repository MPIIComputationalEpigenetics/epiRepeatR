################################################################################
# Perform repeat analysis using the epiRepeatR pipeline
# example call:
# cd /TL/deep/projects/nobackup/fmueller/remoteSync/BroadSVN/eclipse_workspace/repeat_epigenetics/epiRepeatR
# Rscript epiRepeatR.R -a /DEEP_fhgfs/projects/fmueller/repeatEpigenetics/epiRepeatR/analysis/testRuns/tst_run_cl -f testing/data/fileTable_liver.txt --inputBamUnmapped
# testRscript epiRepeatR.R -a /DEEP_fhgfs/projects/fmueller/repeatEpigenetics/epiRepeatR/analysis/testRuns/tst_run_cl2 -f testing/data/fileTable_liver.txt -p 8
################################################################################
suppressPackageStartupMessages(library(epiRepeatR))
library(argparse)

#load the configuration from command line arguments
mainArgs <- loadConfigFromCommandLine()
# print(str(mainArgs))

epiRepeatR:::setUpParallelFromConfig()

#mainArgs <- list(anaDir="/DEEP_fhgfs/projects/fmueller/repeatEpigenetics/epiRepeatR/analysis/testRuns/tst_run_cl2rm",fileTable="testing/data/fileTable_liver.txt")
#epiRepeatR:::setConfigElement("n.processes",8)
res <- do.call("runAnalysis",mainArgs)
