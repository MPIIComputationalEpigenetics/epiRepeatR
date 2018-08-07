library(argparse)

suppressPackageStartupMessages(library(epiRepeatR))
# library(devtools)
# pkg.dir <- "/TL/deep/projects/nobackup/fmueller/remoteSync/BroadSVN/eclipse_workspace/repeat_epigenetics/epiRepeatR"
# pkgName <- file.path(pkg.dir,"epiRepeatR")
# load_all(pkgName)

library(muRtools) #plotting and convenience
library(muLogR) #logging
library(muReportR) #reports

# #testing example
# cmdArgs <- list(
# 	input  = "/DEEP_fhgfs/projects/fmueller/repeatEpigenetics/epiRepeatR/analysis/deepBlood_v04/results/repeatEpigenomeCollection/repeatEpigenomeCollection.rds",
# 	config = "/DEEP_fhgfs/projects/fmueller/repeatEpigenetics/epiRepeatR/analysis/deepBlood_v04/config/config.json",
# 	output = "/DEEP_fhgfs/projects/fmueller/repeatEpigenetics/epiRepeatR/analysis/deepBlood_sandbox/"
# )
# cmdArgs <- list(
# 	input  = "/DEEP_fhgfs/projects/fmueller/repeatEpigenetics/epiRepeatR/analysis/deepBlood_fromGenome_v04mergedInput/results/repeatEpigenomeCollection/repeatEpigenomeCollection.rds",
# 	config = "/DEEP_fhgfs/projects/fmueller/repeatEpigenetics/epiRepeatR/analysis/deepBlood_fromGenome_v04mergedInput/config/config.json",
# 	output = "/DEEP_fhgfs/projects/fmueller/repeatEpigenetics/epiRepeatR/analysis/deepBlood_fromGenome_v04mergedInput/results/exploratoryReport"
# )
# cmdArgs <- list(
# 	input  = "/scratch/users/muellerf/scratch/epiRepeatR/RA_ATAC/RA_repeats_TeTr_ua_v01/results/repeatEpigenomeCollection/repeatEpigenomeCollection.rds",
# 	config = "/scratch/users/muellerf/scratch/epiRepeatR/RA_ATAC/RA_repeats_TeTr_ua_v01/config/config.json",
# 	output = "/scratch/users/muellerf/scratch/epiRepeatR/RA_ATAC/RA_repeats_TeTr_ua_v01/results/exploratoryReport"
# )

ap <- ArgumentParser()
ap$add_argument("-i", "--in", action="store", dest="input", help="Input file (RDS) containing a RepeatEpigenomeCollection object as R dataset.")
ap$add_argument("-o", "--out", dest="output", action="store", help="Output directory")
ap$add_argument("-c", "--config", action="store", help="Config file (json)")
cmdArgs <- ap$parse_args()

logger.cmd.args(cmdArgs)
loadConfig(cmdArgs$config)

if (is.element("debug", names(epiRepeatR:::.config)) && epiRepeatR:::.config$debug){
	saveRDS(cmdArgs, file.path(cmdArgs$output, "cmdargs.rds"))
}

################################################################################
# globals
################################################################################
# replicateColumns <- c("donor","cellType_short","ALL")
replicateColumns <- getConfigElement("annotCols.replicates")
colPal <- colpal.bde

################################################################################
# helpers
################################################################################
plotCorHeatmap <- function(corMat, dendro=FALSE, annotVec=NULL, addCoeff=TRUE, ...){
	paramList <- list(
		corMat,
		Rowv=rev(dendro),
		Colv=dendro,
		symm=FALSE,
		revC=FALSE,
		scale="none",
		symbreaks=TRUE,
		trace="none",
		col=colorpanel(41,"#053061","#f7f7f7","#67001f")
		# margins=c(4,3)
	)
	if (length(annotVec) == ncol(corMat)){
		annotLvls <- as.character(unique(annotVec))
		colorMap <- rep(colPal, length.out=length(annotLvls))
		names(colorMap) <- annotLvls
		colorVec <- colorMap[as.character(annotVec)]
		paramList <- c(paramList, list(ColSideColors=colorVec, RowSideColors=colorVec))
	}
	if (addCoeff){
		paramList <- c(paramList, list(cellnote=round(corMat,2),notecol="grey"))
	}
	paramList <- c(paramList, list(...))
	# print(paramList)
	do.call(heatmap.2, paramList)
}

################################################################################
# main
################################################################################
theme_set(theme_bw())

logger.start("Preparing data")
	rec <- readRDS(cmdArgs$input)
	resDir <- cmdArgs$output

	if (!dir.exists(resDir)){
		dir.create(resDir, recursive=TRUE)
		logger.info(c("Created output directory:",resDir))
	}

	markLvls <- getMarks(rec)
	markLvls.named <- markLvls
	names(markLvls.named) <- markLvls

	logger.start("Filtering")
		if (is.element("DNAmeth", markLvls)){
			rec <- filterRepRefMeth(
				rec,
				minReads=getConfigElement("plotRepTree.meth.minReads"),
				minCpGs=getConfigElement("plotRepTree.meth.minCpGs"),
				minCpGcov=getConfigElement("meth.minCpGcov")
			)
		}
		if (length(setdiff(markLvls, "DNAmeth")) > 0){
			rec <- filterRepRefChip(rec, minReads=getConfigElement("plotRepTree.meth.minReads"))
			rec <- filterRepRefAcc(rec, minReads=getConfigElement("plotRepTree.meth.minReads"))
		}
		numReps <- length(epiRepeatR:::getRepeatIds(getRepRef(rec)))
		logger.info(c(numReps, "repeats retained after filtering"))
	logger.completed()

	aa <- getAnnot(rec)
	aa.all <- data.frame(aa, ALL="all")

	sampleGroups <- getSampleGroups(aa, addAll=TRUE)
	sampleGroupNames <- names(sampleGroups)
	names(sampleGroupNames) <- normalize.str(sampleGroupNames, resolve.camel=TRUE, return.camel=TRUE)

	report <- createReport(
		file.path(resDir, "repeatExploration.html"),
		"Repeat Exploratory Analysis",
		page.title="epiRepeatR Report",
		init.configuration=TRUE
	)
logger.completed()

logger.start("Retrieving score matrices")
	scoreMatList <- lapply(markLvls, FUN=function(mn){
		res <- getRepeatScores(rec, mn, dropEmptySamples=TRUE)
		res[res==Inf] <- max(res[res!=Inf], na.rm=TRUE)
		res[res==-Inf] <- min(res[res!=-Inf], na.rm=TRUE)
		return(res)
	})
	names(scoreMatList) <- markLvls
	samples4marks <- lapply(markLvls, FUN=function(mn){
		sns <- intersect(colnames(scoreMatList[[mn]]), rownames(aa.all))
		if (length(sns) < nrow(aa.all)){
			naSamples <- setdiff(rownames(aa.all), sns)
			logger.warning(c("The following samples do not have observations (",mn,"):", paste(naSamples, collapse=",")))
		}
		return(sns)
	})
	names(samples4marks) <- markLvls
logger.completed()

logger.start("Dimension reduction plots")
	dimRedMethods <- c(
		"pca"="PCA",
		"mdseuc"="MDS (euc)", "mdsman"="MDS (man)"
		# "tsneeuc"="tSNE (euc)", "tsneman"="tSNE (man)"
	)
	plotList <- list()
	for (mn in markLvls){
		logger.start(c("Mark:", mn))
			aaCur <- aa.all[samples4marks[[mn]],]
			scoreMatT <- t(scoreMatList[[mn]])

			coordList <- list()
			coordList[["pca"]]     <- getDimRedCoords.pca(scoreMatT)
			coordList[["mdseuc"]]  <- getDimRedCoords.mds(scoreMatT, distMethod="euclidean")
			coordList[["mdsman"]]  <- getDimRedCoords.mds(scoreMatT, distMethod="manhattan")
			# coordList[["tsneeuc"]] <- getDimRedCoords.tsne(scoreMatT, distMethod="euclidean")
			# coordList[["tsneman"]] <- getDimRedCoords.tsne(scoreMatT, distMethod="manhattan")

			for (i in 1:length(sampleGroupNames)){
				for (j in 1:length(sampleGroupNames)){
					ppl <- lapply(names(dimRedMethods), FUN=function(drm){
						res <- list(
							plot = getDimRedPlot(
								coordList[[drm]],
								annot=aaCur,
								colorCol=sampleGroupNames[i],
								shapeCol=sampleGroupNames[j],
								colScheme=colPal,
								addLabels=FALSE,
								addDensity=FALSE,
								annot.text=NULL
							),
							fileName = paste("dimRed", drm, mn, names(sampleGroupNames)[i], names(sampleGroupNames)[j], sep="_")
						)
						return(res)
					})
					plotList <- c(plotList, ppl)
				}
				#add a plot containing the labels instead of shape
				ppl <- lapply(names(dimRedMethods), FUN=function(drm){
					res <- list(
						plot = getDimRedPlot(
							coordList[[drm]],
							annot=aaCur,
							colorCol=sampleGroupNames[i],
							shapeCol="ALL",
							colScheme=colPal,
							addLabels=TRUE,
							addDensity=FALSE,
							annot.text=NULL
						),
						fileName = paste("dimRed", drm, mn, names(sampleGroupNames)[i], "labels", sep="_")
					)
					return(res)
				})
				plotList <- c(plotList, ppl)
			}
		logger.completed()
	}

	stext <- "The following plot shows dimension reduction for the epigenetic marks in repeats:"
	report <- addReportSection(report, "Dimension Reduction", stext)

	logger.start("Plotting")
		figPlots <- list()
		for (ple in plotList){
			rplot <- createReportGgPlot(ple$plot, ple$fileName, report, create.pdf=TRUE, width=7, height=7, high.png=200)
			rplot <- off(rplot, handle.errors=TRUE)
			figPlots <- c(figPlots,list(rplot))
		}
	logger.completed()

	settingNames <- list(
		"Method"=dimRedMethods,
		"Mark"=markLvls.named,
		"Color by"=sampleGroupNames,
		"Shape by"= c(sampleGroupNames, "labels"="label")
	)
	#default figure
	defIndex <- grep(paste("dimRed", "pca", markLvls[1], "all", "all", sep="_"), sapply(figPlots, FUN=function(x){x@fname}))[1]
	if (length(defIndex) < 1) defIndex <- 1
	desc <- "Dimension reduction plots for epigenetic marks in repeats"
	report	<- addReportFigure(report, desc, figPlots, settingNames, selected.image=defIndex)
logger.completed()

logger.start("Sample Correlation Matrix")
	figPlots <- list()
	for (mn in markLvls){
		logger.start(c("Mark:", mn))
			aaCur <- aa.all[samples4marks[[mn]],]
			
			corMat <- cor(scoreMatList[[mn]], use="pairwise.complete.obs")
			dd <- as.dist(1 - corMat)
			cres <- hclust(dd, method="ward.D2")
			cresDend <- as.dendrogram(cres)

			for (i in 1:length(sampleGroupNames)){
				annotVec <- aaCur[,sampleGroupNames[i]]

				fName <- paste("corHm", mn, names(sampleGroupNames)[i], sep="_")

				rplot <- createReportPlot(fName, report, create.pdf=TRUE, width=7, height=7, high.png=200)
				plotCorHeatmap(corMat, dendro=cresDend, annotVec=annotVec)
				off(rplot)
				figPlots <- c(figPlots, list(rplot))
			}
		logger.completed()
	}

	stext <- "The following plot shows sample correlation coefficients for the epigenetic marks in repeats:"
	report <- addReportSection(report, "Sample Correlation", stext)
	settingNames <- list(
		"Mark"=markLvls.named,
		"Color by"=sampleGroupNames
	)
	#default figure
	defIndex <- grep(paste("corHm", markLvls[1], "all", sep="_"), sapply(figPlots, FUN=function(x){x@fname}))[1]
	if (length(defIndex) < 1) defIndex <- 1
	desc <- "Correlation heatmap for epigenetic marks in repeats. Pearson correlation is shown."
	report	<- addReportFigure(report, desc, figPlots, settingNames, selected.image=defIndex)
logger.completed()

if (length(replicateColumns) > 0){
	logger.start("Replicate Scatterplots")
		repPairList <- list()
		for (rc in replicateColumns){
			annotF <- factor(aa.all[,rc])
			annotT <- table(annotF)
			annotT <- annotT[annotT>1]
			repIndexPairs <- list()
			lvls <- c()
			if (length(annotT) > 0){
				lvls <- names(annotT)
				for (curLvl in lvls){
					lvlInds <- which(annotF==curLvl)
					lvlPairs <- combn(lvlInds, 2)
					rl <- list(
						indices=lvlInds,
						indexPairs=lvlPairs,
						numPairs=ncol(lvlPairs),
						annotColumn=rc,
						level=curLvl
					)
					repPairList <- c(repPairList, list(rl))
				}
			}
		}

		plotList <- list()
		for (mn in markLvls){
			logger.start(c("Mark:", mn))
				scoreMat <- getRepeatScores(rec, mn)
				scoreMat[scoreMat==Inf] <- max(scoreMat[scoreMat!=Inf], na.rm=TRUE)
				scoreMat[scoreMat==-Inf] <- min(scoreMat[scoreMat!=-Inf], na.rm=TRUE)
				for (i in 1:length(repPairList)){
					lvlPairs <- repPairList[[i]]$indexPairs
					for (j in 1:ncol(lvlPairs)){
						df2p <- data.frame(scoreMat[,lvlPairs[,j]])
						cc <- cor(df2p[,1],df2p[,2],use="pairwise.complete.obs")
						#if one or both of the samples is not present for the current mark
						#cc will be NA and thus no scatterplot can be generated
						pp <- ggMsgPlot("N/A")
						if (!is.na(cc)){
							pp <- create.densityScatter(df2p)
						}
						plotList <- c(plotList, list(list(
							plot=pp,
							fileName=paste("repScatter", mn, paste0("c",i,"p",j), sep="_")
						)))
					}
				}
			logger.completed()
		}

		stext <- "The following plot shows the agreement between replicates according to annotation:"
		report <- addReportSection(report, "Replicate Agreement", stext)

		logger.start("Plotting")
			figPlots <- list()
			for (ple in plotList){
				rplot <- createReportGgPlot(ple$plot, ple$fileName, report, create.pdf=TRUE, width=7, height=7, high.png=200)
				rplot <- off(rplot, handle.errors=TRUE)
				figPlots <- c(figPlots,list(rplot))
			}
		logger.completed()

		#construct the names of the comparisons
		sns <- rownames(aa.all)
		compVec <- c()
		compVecNames <- c()
		for (i in 1:length(repPairList)){
			cvnsCur <- paste0("c",i, "p", 1:repPairList[[i]]$numPairs)
			compVecNames <- c(compVecNames, cvnsCur)

			sampleNamePairs <- apply(repPairList[[i]]$indexPairs, 2, FUN=function(pp){
				paste(sns[pp], collapse=" vs. ")
			})
			cvsCur <- paste0(repPairList[[i]]$annotColumn, " > ", repPairList[[i]]$level, " > ", sampleNamePairs)
			compVec <- c(compVec, cvsCur)
		}
		names(compVec) <- compVecNames

		settingNames <- list(
			"Mark"=markLvls.named,
			"Comparison"=compVec
		)
		desc <- "Scatterplot of replicate comparisons"
		report	<- addReportFigure(report, desc, figPlots, settingNames)
	logger.completed()
}

off(report)

################################################################################
# sandbox
################################################################################
if (FALSE){

aa$donor <- factor(paste0("HX0",aa$donor))
repMeth <- getRepeatScores(rec, "DNAmeth")

#dimension reduction
ppl <- plotAllDimRed(
	t(repMeth),
	file.path(resDir, "repMeth_dimRed"), fn.suffix="",
	annot=aa,
	colorCol="cellType_short",
	shapeCol="donor",
	colScheme=colpal.bde
)

plotReplicateScatters <- function(scoreMat, repAnnot, ...){
	annotF <- factor(repAnnot)
	annotT <- table(annotF)
	annotT <- annotT[annotT>1]
	res <- list()
	if (length(annotT) > 0){
		for (curLvl in names(annotT)){
			lvlInds <- which(annotF==curLvl)
			lvlPairs <- combn(lvlInds, 2)
			for (i in 1:ncol(lvlPairs)){
				df2p <- data.frame(scoreMat[,lvlPairs[,i]])
				cc <- cor(df2p[,1],df2p[,2],use="pairwise.complete.obs")
				pp <- create.densityScatter(df2p, ...)
				res <- c(res, list(list(
					plot=pp,
					correlation=cc,
					level=curLvl,
					samples=colnames(scoreMat)[lvlPairs[,i]]
				)))
			}
		}
	}
	return(res)
}
ppl <- plotReplicateScatters(repMeth, aa$cellType_short, sparse.points=0.01, add.text.cor=TRUE)
if (length(ppl) > 0){
	for (i in 1:length(ppl)){
		fName <- paste0("repScatter_", "cellType", "_", ppl[[i]]$level, "_", ppl[[i]]$samples[1], "_vs_", ppl[[i]]$samples[2], ".pdf")
		fName <- file.path(resDir, fName)
		ggsave(fName, ppl[[i]]$plot, width=10, height=10)
	}
}

}
