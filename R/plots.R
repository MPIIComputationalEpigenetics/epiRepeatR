#' plotHeatmap
#'
#' Given a metrix of values, plots a corresponding heatmap into the assigned plot space. Optionally adds text containing the values to the boxes
#'
#' @param X				a matrix of values to be plotted
#' @param xmin			bounding coordinates for the heatmap
#' @param xmax 			bounding coordinates for the heatmap
#' @param ymin			bounding coordinates for the heatmap
#' @param ymax			bounding coordinates for the heatmap
#' @param colorGradient	a vector of colors to be used for the heatmap (low values to high values)
#' @param zlim			limits to be applied to match the colors to values in the matrix
#' @param addValueTxt	flag indicating whether values should be added as text to the heatmap boxes
#' @param roundDigits.valTxt	if text is added, to how many digits should the values be rounded
#' @param cex.valTxt	text size for the value text in the boxes
#' @param col.valTxt	text color for the value text in the boxes
#' @return an object (S3) containing coordinates of gridlines and box centers of the added heatmap
#' @author Fabian Mueller
#' @noRd
plotHeatmap <- function(X, xmin=0, xmax=1, ymin=0, ymax=1, colorGradient=colorpanel(100,"#EDF8B1","#41B6C4","#081D58"), zlim=c(min(X,na.rm=TRUE),max(X,na.rm=TRUE)), addValueTxt=FALSE, roundDigits.valTxt=2, cex.valTxt=1.0, col.valTxt="black"){
	gridlinePos.x <- seq(xmin, xmax, length.out=ncol(X)+1) #note that the matrix will be transposed for plotting later. Hence switch rows and columns
	gridlinePos.y <- seq(ymin, ymax, length.out=nrow(X)+1)

	gridMidPoints.x <- gridlinePos.x[1:(length(gridlinePos.x)-1)] + (gridlinePos.x[2:length(gridlinePos.x)] - gridlinePos.x[1:(length(gridlinePos.x)-1)]) / 2
	gridMidPoints.y <- gridlinePos.y[1:(length(gridlinePos.y)-1)] + (gridlinePos.y[2:length(gridlinePos.y)] - gridlinePos.y[1:(length(gridlinePos.y)-1)]) / 2
	
	matrixMidPoints.x <- matrix(rep(gridMidPoints.x, nrow(X)), ncol=ncol(X), byrow=TRUE)
	matrixMidPoints.y <- matrix(rep(gridMidPoints.y, ncol(X)), ncol=ncol(X), byrow=FALSE)
	X.trunc <- t(X)
	X.trunc[X.trunc < zlim[1]] <- zlim[1]
	X.trunc[X.trunc > zlim[2]] <- zlim[2]
	image(
		gridlinePos.x, gridlinePos.y,
		X.trunc,
		zlim=zlim,
		axes=FALSE, col=colorGradient, add=TRUE
	)
	if (addValueTxt){
		valTxt <- sprintf(paste0("%.",roundDigits.valTxt,"f"), round(X, roundDigits.valTxt))
		text(as.vector(matrixMidPoints.x), as.vector(matrixMidPoints.y), labels=valTxt, cex=cex.valTxt, col=col.valTxt)
	}
	res <- list(
		gridlinePos.x=gridlinePos.x,
		gridlinePos.y=gridlinePos.y,
		gridMidPoints.x=gridMidPoints.x,
		gridMidPoints.y=gridMidPoints.y
	)
	class(res) <- "plotHeatmapResult"
	return(res)
}

#' repPlot_groupSummary
#'
#' Plots a repeat sequence group summary tree
#'
#' @param repRef		repeat reference. Object of type \code{\linkS4class{RepeatReference}}.
#' @param scores		a matrix of scores to be plotted. should have one row for each repeat element in the reference and a column for each sample/experiment
#' @param grpInfo 		list of assignments of samples/experiments to groups. Obtained by the \code{\link{getSampleGroups}} function
#' @param colorGradient	a vector of colors to be used for the heatmap (low values to high values). Alternatively, a list of vectors containing one element for each group.
#' @param zlim			vector of length 2 containing the limits to be applied to match the colors to values in the matrix. Alternatively, a list of vectors containing one element for each group.
#' @param groupColors	vector of colors to be used for the sample groups
#' @param txt.cex		base text size to be used in the plot
#' @param leafColors	colors of the leaf nodes/repeat elements in the same order as in \code{getRepeatIds(repRef)}. set to \code{NULL} (default) to disable custom leaf color
#' @param dendroMethod  method for plotting the repeat subfamily dendrogram. See \code{RepeatTree} class for possible values.
#' @return nothing of particular interest
#'
#' @details
#' Best see for yourself what the result looks like. The resulting plot contains a dendrogram of repetitive element, groupwise heatmaps
#' of supplied values, group averages of supplied values and per-group boxplots for each repetitive element.
#'
#' @author Fabian Mueller
#' @noRd
repPlot_groupSummary <- function(
		repRef,
		scores,
		grpInfo,
		colorGradient=colorpanel(100,"#EDF8B1","#41B6C4","#081D58"),
		zlim=c(min(scores, na.rm=TRUE), max(scores, na.rm=TRUE)),
		groupColors=rep(c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#2166AC","#B2182B","#00441B","#40004B","#053061"), length.out=length(grpInfo)),
		txt.cex=0.5,
		leafColors=NULL,
		addBoxplots=FALSE,
		dendroMethod="repeatFamily"){

	suppressPackageStartupMessages(require(diagram))

	repTree <- RepeatTree(repRef, method=dendroMethod)
	# repTree <- RepeatTree(repRef, method="hierClust")
	orderedLabels       <- getDendrogramMembers(repTree, rev=FALSE)
	orderedLabels4plot  <- getDendrogramMembers(repTree, rev=TRUE)
	scores.ordered4plot <- scores[orderedLabels4plot,,drop=FALSE]
	
	#setting the coordinates for the elements
	Ngrps <- length(grpInfo)
	Nseqs <- nrow(getRepeatInfo(repRef))

	# convert arguments that can differ for each group to group lists
	if (!is.list(colorGradient))  colorGradient <- list(colorGradient)
	if (length(colorGradient)==1) colorGradient <- rep(colorGradient, Ngrps)
	if (length(colorGradient) != Ngrps) stop(paste0("Invalid argument: colorGradient. Must be of length 1 or N_groups"))
	if (!is.list(zlim))  zlim <- list(zlim)
	if (length(zlim)==1) zlim <- rep(zlim, Ngrps)
	if (length(zlim) != Ngrps) stop(paste0("Invalid argument: zlim. Must be of length 1 or N_groups"))
	#aggregate zlim for combining all groups
	zlim.c <- c(min(unlist(zlim)), max(unlist(zlim)))		


	coord.tree <- c(xmin=0.0, xmax=0.15, ymin=0.05, ymax=0.95)
	grpHeat.start.x <- 0.18
	grpHeat.ymin <- coord.tree["ymin"]
	grpHeat.ymax <- coord.tree["ymax"]

	txt.cex.base <- txt.cex
	spacer.w <- 0.01 # custom width for a spacer
	heatBox.w <- 0.015
	heatBox.h <- (coord.tree["ymax"] - coord.tree["ymin"]) / Nseqs #height of a heatmap box

	#for each group get the coordinates for their heatmap
	coord.heat.grps <- matrix(NA,nrow=4,ncol=Ngrps)
	rownames(coord.heat.grps) <- c("xmin","xmax","ymin","ymax")
	prevEnd.x <- grpHeat.start.x
	for (i in 1:Ngrps){
		start.x <- prevEnd.x + spacer.w
		add.x <- heatBox.w * length(grpInfo[[i]])
		coord.heat.grps[,i] <- c(start.x, start.x + add.x, grpHeat.ymin, grpHeat.ymax)
		prevEnd.x <- start.x + add.x
	}
	#group summary heatmap coordinates
	coord.heat.grp.sum <- matrix(NA,nrow=4,ncol=Ngrps)
	rownames(coord.heat.grp.sum) <- c("xmin","xmax","ymin","ymax")
	coord.heat.grp.sum["xmin",] <- max(coord.heat.grps["xmax",]) + spacer.w + heatBox.w * 0:(Ngrps-1)
	coord.heat.grp.sum["xmax",] <- coord.heat.grp.sum["xmin",] + heatBox.w
	coord.heat.grp.sum["ymin",] <- grpHeat.ymin
	coord.heat.grp.sum["ymax",] <- grpHeat.ymax

	if (addBoxplots){
		bp.w <- 4 * heatBox.w #width of the boxplot
		bp.box.h <- ((coord.tree["ymax"] - coord.tree["ymin"]) / Nseqs) * 0.90
		coord.bp <- c(
			xmin=max(coord.heat.grp.sum["xmax",])+spacer.w,
			xmax=max(coord.heat.grp.sum["xmax",])+spacer.w+bp.w,
			ymin=unname(grpHeat.ymin),
			ymax=unname(grpHeat.ymax)
		)
	}

	#set the color values for the repeat tree leafs if required
	if (!is.null(leafColors)){
		leafColors.ordered <- leafColors
		names(leafColors.ordered) <- getRepeatIds(repRef)
		leafColors.ordered <- leafColors.ordered[orderedLabels]
		repTree <- setDendrogramMemberAttr(repTree, "leafBgColor", leafColors.ordered)
	}

	#start plotting
	openplotmat()

	#plot the repeat tree
	addToPlot(repTree, coord.tree["xmin"], coord.tree["xmax"], coord.tree["ymin"], coord.tree["ymax"])

	#for each group plot a heatmap
	for (i in 1:Ngrps){
		scoreMat <- scores.ordered4plot[,grpInfo[[i]], drop=FALSE]
		heatRes <- plotHeatmap(
			scoreMat,
			coord.heat.grps["xmin",i], coord.heat.grps["xmax",i], coord.heat.grps["ymin",i], coord.heat.grps["ymax",i],
			colorGradient=colorGradient[[i]],
			zlim=zlim[[i]],
			addValueTxt=TRUE,
			roundDigits.valTxt=2,
			cex.valTxt=txt.cex.base,
			col.valTxt="dark grey"
		)
		#plot the column headers
		text(heatRes$gridMidPoints.x, coord.heat.grps["ymax",i], grpInfo[[i]], srt=45,adj=c(0,0), cex=txt.cex.base, col=groupColors[i])

		#plot the group summary heatmap
		summaryColNames <- paste0(names(grpInfo)[i],"_mean")
		scoreMatMean <- matrix(rowMeans(scoreMat, na.rm=TRUE))
		rownames(scoreMatMean) <- rownames(scoreMat)
		colnames(scoreMatMean) <- summaryColNames
		heatResMean <- plotHeatmap(
			scoreMatMean,
			coord.heat.grp.sum["xmin",i], coord.heat.grp.sum["xmax",i], coord.heat.grp.sum["ymin",i], coord.heat.grp.sum["ymax",i],
			colorGradient=colorGradient[[i]],
			zlim=zlim[[i]],
			addValueTxt=TRUE,
			roundDigits.valTxt=2,
			cex.valTxt=txt.cex.base,
			col.valTxt="dark grey"
		)
		text(heatResMean$gridMidPoints.x, coord.heat.grp.sum["ymax",i], summaryColNames, srt=45, adj=c(0,0), cex=txt.cex.base, col=groupColors[i])

		#plot groupwise boxplots
		if (addBoxplots){
			#linearly scale values to fit in the interval for the boxplots
			bp.data <- lapply(1:nrow(scoreMat), FUN=function(i){
				x <- scoreMat[i,]
				return((x-zlim.c[1]) / (zlim.c[2]-zlim.c[1]) * (coord.bp["xmax"]-coord.bp["xmin"]) + coord.bp["xmin"])
			})
			boxplot(bp.data,horizontal=T,add=T,at=heatRes$gridMidPoints.y, border=groupColors[i],boxwex=bp.box.h,axes=FALSE)
		}
	}

	if (addBoxplots){
		#axis for boxplots
		tickmarks <- seq(coord.bp["xmin"],coord.bp["xmax"],length.out=3)
		axis.txt <- sprintf(paste0("%.",2,"f"), round(c(zlim.c[1],mean(zlim.c),zlim.c[2]), 2))
		axis(3,at=tickmarks,pos=coord.bp["ymax"],labels=axis.txt, cex.axis=txt.cex)
		for (tm in tickmarks){
			lines(x=c(tm,tm), y=coord.bp[c("ymin","ymax")],col="grey",lty="dashed")
		}
	}
}


#' createRepPlot_groupSummaryTrees_meth
#'
#' Plots a repeat sequence group summary tree for methylation values for each grouping information entry
#'
#' @param .obj	            \code{\linkS4class{RepeatEpigenomeCollection}} object
#' @param plotDir			Output directory where the plots are saved to
#' @param dendroMethod      method for plotting the repeat subfamily dendrogram. See \code{RepeatTree} class for possible values.
#' @param minReads			filtering parameter: specify the minimum number of reads that must match to a given repeat element in order for the repeat to be added to the plot
#' @param minCpGs			filtering parameter: specify the minimum number of CpG that must be contained in a given repeat element in order for the repeat to be added to the plot
#' @param minCpGcov			filtering parameter: specify the minimum number of reads that must cover a given CpG in order to be considered a valid methylation measurement
#' @return nothing of particular interest
#'
#' @details
#' Which repeats are used is regulated by the filtering parameters.
#' One PDF file for each comparison information is created in the output directory
#'
#' @author Fabian Mueller
#' @noRd
createRepPlot_groupSummaryTrees_meth <- function(
		.obj,
		plotDir,
		dendroMethod="repeatFamily",
		minReads=getConfigElement("plotRepTree.meth.minReads"),
		minCpGs=getConfigElement("plotRepTree.meth.minCpGs"),
		minCpGcov=getConfigElement("meth.minCpGcov")){
	suppressPackageStartupMessages(require(gplots))
	suppressPackageStartupMessages(require(diagram))

	rec <- filterRepRefMeth(.obj, minReads=minReads, minCpGs=minCpGs, minCpGcov=minCpGcov)
	repRef <- getRepRef(rec)
	if (length(getRepeatIds(repRef)) < 1) logger.error("No repeat sequence retained after filtering")

	methScores <- getRepeatScores(rec, "DNAmeth", dropEmptySamples=TRUE, minCpGcov=minCpGcov)
	covgMat    <- getRepeatCovg(rec, "DNAmeth", dropEmptySamples=TRUE)
	if (ncol(covgMat) < 1){
		logger.info("No read coverage found. --> using the number of repeat instances as a proxy.")
		covgMat <- getRepeatCovg(rec, "DNAmeth", dropEmptySamples=TRUE, type="numInstances")
	}
	sampleNames <- colnames(methScores)
	if (any(sampleNames != colnames(covgMat))){
		logger.error("Nonmatching sample names for meth scores and coverage")
	}

	annot <- getAnnot(rec)
	rownames(annot) <- getSamples(rec)
	annot <- annot[sampleNames,]
	sampleGroups <- getSampleGroups(annot, addAll=TRUE)

	#create a plot for each grouping info
	for (i in 1:length(sampleGroups)){
		ggn <- names(sampleGroups)[i]
		ggs <- sampleGroups[[i]]

		covgMat.sub <- covgMat[,unlist(ggs)]
		covgMat.sub.rel <-  apply(covgMat.sub, 2, FUN=function(x){x/sum(x)}) # relative coverage
		leafColScore <- rowMeans(covgMat.sub.rel, na.rm=TRUE)
		#order the color values
		leafColors <- colorize.value(leafColScore, colscheme=colorpanel(100,"white","black"))
		fn <- file.path(plotDir, paste0("repeatTree_groupSummary_",ggn,".pdf"))
		# load_all(file.path(pkg.dir,"epiRepeatR"))
		pdf(fn, width=20, height=100)
			repPlot_groupSummary(repRef, methScores, ggs, zlim=c(0,1), leafColors=leafColors, addBoxplots=TRUE, dendroMethod=dendroMethod)
		dev.off()
	}
	res <- list(
		status="success",
		callParams=list(
			plotDir=plotDir,
			minReads=minReads,
			minCpGs=minCpGs
		)
	)
	invisible(res)
}


#' createRepPlot_markTree
#'
#' Plots a repeat sequence summary tree for all marks in the dataset
#'
#' @param .obj	            \code{\linkS4class{RepeatEpigenomeCollection}} object
#' @param plotDir			Output directory where the plots are saved to
#' @param dendroMethod      method for plotting the repeat subfamily dendrogram. See \code{\linkS4class{RepeatTree-class}} class for possible values.
#' @param leafColorMethod   method for coloring the leafs of the dendrogram. Options are \code{"coverage"} for read coverage (default) and \code{"abundance"} for genomic abundance
#' @param normChipMethod    method for normalizing ChIP enrichment scores per mark. See \code{\link{normalizeMatrix}} for possible values. 
#' @param minReads			filtering parameter: specify the minimum number of reads that must match to a given repeat element in order for the repeat to be added to the plot
#' @param minCpGs			filtering parameter: specify the minimum number of CpG that must be contained in a given repeat element in order for the repeat to be added to the plot
#' @param minCpGcov			filtering parameter: specify the minimum number of reads that must cover a given CpG in order to be considered a valid methylation measurement
#' @return nothing of particular interest
#'
#' @details
#' The reference repeats are obtained via \code{RepeatReference()}. Which repeats are used is regulated by the filtering parameters
#' One PDF file for each comparison information is created in the output directory
#'
#' @author Fabian Mueller
#' @noRd
createRepPlot_markTree <- function(
		.obj,
		plotDir,
		dendroMethod="repeatFamily",
		leafColorMethod=getConfigElement("plotRepTree.leafColorMethod"),
		normChipMethod="none",
		minReads=getConfigElement("plotRepTree.meth.minReads"),
		minCpGs=getConfigElement("plotRepTree.meth.minCpGs"),
		minCpGcov=getConfigElement("meth.minCpGcov")){
	suppressPackageStartupMessages(require(gplots))
	suppressPackageStartupMessages(require(diagram))

	rec <- filterRepRefMeth(.obj, minReads=minReads, minCpGs=minCpGs, minCpGcov=minCpGcov)
	rec <- filterRepRefChip(rec, minReads=minReads)
	rec <- filterRepRefAcc(rec, minReads=minReads)
	repRef <- getRepRef(rec)
	if (length(getRepeatIds(repRef)) < 1) logger.error("No repeat sequence retained after filtering")

	sampleNames <- getSamples(rec)
	markLvls <- getMarks(rec)

	sampleGroups <- getSampleGroups(getAnnot(rec), addAll=TRUE)

	sampleNames.exp <- rep(sampleNames, length(markLvls))
	markNames.exp   <- rep(markLvls, each=length(sampleNames))
	sampleNames.full <- paste(sampleNames.exp, markNames.exp, sep="_")

	scoreMat <- do.call("cbind", lapply(markLvls, FUN=function(mn){
		rr <- getRepeatScores(rec, mn, minCpGcov=minCpGcov)
		mt <- inferMarkTypes(mn)
		if (is.element(mt, c("ChIPseq", "Acc"))){
			rr <- normalizeMatrix(rr, method=normChipMethod)
		}
		return(rr)
	}))
	colnames(scoreMat) <- sampleNames.full
	hasScore <- apply(scoreMat,2,FUN=function(x){!all(is.na(x))})

	covgMat <- do.call("cbind", lapply(markLvls, FUN=function(mn){
		cm <- getRepeatCovg(rec, mn)
		if (mn == "DNAmeth" && all(is.na(cm))){
			logger.info("No read coverage found. --> using the number of repeat instances as a proxy.")
			cm <- getRepeatCovg(rec, mn, type="numInstances")
		}
		return(cm)
	}))
	colnames(covgMat) <- sampleNames.full

	leafColScore <- NULL
	if (leafColorMethod=="abundance"){
		repInfo <- getRepeatInfo(repRef)
		leafColScore <- repInfo[,"percCovg"]
	} else {
		covgMat.rel <- apply(covgMat, 2, FUN=function(x){x/sum(x)}) # relative coverage
		#TODO: add the option to color leaves by genomic abundance (if that annotation is in the REC .obj)
		leafColScore <- rowMeans(covgMat.rel, na.rm=TRUE)
	}
	leafColors <- colorize.value(leafColScore, colscheme=colorpanel(100,"white","black"))
	leafColors[is.na(leafColScore)] <- "#bd0026" #"red"

	markGroups <- lapply(markLvls, FUN=function(mn){
		sampleNames.full[hasScore & markNames.exp==mn]
	})
	names(markGroups) <- markLvls

	colGrads <- rep(list(colorpanel(100,"#8C510A","#F5F5F5","#01665E")), length(markLvls))
	zlims <- rep(list(c(-3, 3)), length(markLvls)) #chip-seq abs(log2FC) limit to range [-3,3]
	indMeth <- match("DNAmeth", markLvls)
	if (!is.na(indMeth)){
		colGrads[[indMeth]] <- colorpanel(100,"#EDF8B1","#41B6C4","#081D58")
		zlims[[indMeth]] <- c(0,1)
	}
	indAcc <- match("Acc", inferMarkTypes(markLvls))
	if (!is.na(indAcc) && length(indAcc)){
		for (i in indAcc){
			zlims[[i]] <- c(-5,5)
		}
	}

	fn <- file.path(plotDir, paste0("repeatTree_markSummary",".pdf"))
	pdf(fn, width=20, height=100)
		repPlot_groupSummary(repRef, scoreMat, markGroups, colorGradient=colGrads, zlim=zlims, leafColors=leafColors, addBoxplots=FALSE, dendroMethod=dendroMethod)
	dev.off()

	if (length(sampleGroups) > 0){
		for (i in 1:length(sampleGroups)){
			ggn <- names(sampleGroups)[i]
			ggs <- sampleGroups[[i]]
			ggs.sn <- lapply(ggs, FUN=function(gg){
				sampleNames.full[sampleNames.exp %in% gg]
			})
			names(ggs.sn) <- names(ggs)

			scoreMat.g <- do.call("cbind", lapply(markLvls, FUN=function(mn){
				mm <- do.call("cbind", lapply(names(ggs), FUN=function(gn){
					curSampleNames <- intersect(ggs.sn[[gn]], markGroups[[mn]])
					rowMeans(scoreMat[, curSampleNames, drop=FALSE])
				}))
				colnames(mm) <- paste(names(ggs), mn, sep="_")
				return(mm)
			}))
			ggs.mark <- lapply(markLvls, FUN=function(mn){
				paste(names(ggs), mn, sep="_")
			})
			names(ggs.mark) <- markLvls

			if (leafColorMethod!="abundance"){
				covgMat.sub <- covgMat[,unlist(ggs.sn)]
				covgMat.sub.rel <-  apply(covgMat.sub, 2, FUN=function(x){x/sum(x)}) # relative coverage
				leafColScore <- rowMeans(covgMat.sub.rel, na.rm=TRUE)
				#order the color values
				leafColors <- colorize.value(leafColScore, colscheme=colorpanel(100,"white","black"))
				leafColors[is.na(leafColScore)] <- "#bd0026" #"red"
			}

			fn <- file.path(plotDir, paste0("repeatTree_markSummary_",ggn,".pdf"))
			pdf(fn, width=20, height=100)
				repPlot_groupSummary(repRef, scoreMat.g, ggs.mark, colorGradient=colGrads, zlim=zlims, leafColors=leafColors, addBoxplots=FALSE, dendroMethod=dendroMethod)
			dev.off()
		}
	}

	res <- list(
		status="success",
		callParams=list(
			plotDir=plotDir,
			minReads=minReads,
			minCpGs=minCpGs
		)
	)
}

################################################################################
# Differential plots
################################################################################

#' @param repRef		repeat reference. Object of type \code{\linkS4class{RepeatReference}}.
#' @param compInfo      list of comparison info as returned by \code{\link{getComparisonInfo,RepeatEpigenomeCollection-method}}
#' @param diffScores	list of differential score matrices. Obtained by the \code{\link{getRepeatScoresDiff,RepeatEpigenomeCollection-method}} function
#' @param sampleScores	list of sample score matrices. Can for instance be obtained by the \code{\link{getRepeatScores,RepeatEpigenomeCollection-method}} function
#' @param colorGradient.groupScore	a vector of colors to be used for the group score heatmap (low values to high values). Alternatively, a list of vectors containing one element for each comparison/mark.
#' @param zlim.groupScore			vector of length 2 containing the limits to be applied to match the colors to values in the group score matrix. Alternatively, a list of vectors containing one element for each comparison/mark.
#' @param colorGradient.diff	a vector of colors to be used for the differential score heatmap (low values to high values). Alternatively, a list of vectors containing one element for each comparison/mark.
#' @param zlim.diff			vector of length 2 containing the limits to be applied to match the colors to values in the differential score matrix. Alternatively, a list of vectors containing one element for each comparison/mark.
#' @param colorGradient.score	a vector of colors to be used for the sample repeat score heatmap (low values to high values). Alternatively, a list of vectors containing one element for each comparison/mark.
#' @param zlim.score			vector of length 2 containing the limits to be applied to match the colors to values in the sample repeat score matrix. Alternatively, a list of vectors containing one element for each comparison/mark.
#' @param leafColors	colors of the leaf nodes/repeat elements in the same order as in \code{getRepeatIds(repRef)}. set to \code{NULL} (default) to disable custom leaf color
#' @param dendroMethod  method for plotting the repeat subfamily dendrogram. See \code{RepeatTree} class for possible values.
repPlot_differential <- function(
		repRef,
		compInfo,
		diffScores,
		sampleScores=list(),
		colorGradient.groupScore=colorpanel(100,"#EDF8B1","#41B6C4","#081D58"),
		zlim.groupScore=NULL,
		colorGradient.diff=colorpanel(100,"#EDF8B1","#41B6C4","#081D58"),
		zlim.diff=NULL,
		colorGradient.score=colorpanel(100,"#EDF8B1","#41B6C4","#081D58"),
		zlim.score=NULL,
		leafColors=NULL,
		dendroMethod="repeatFamily"){

	suppressPackageStartupMessages(require(ComplexHeatmap))
	suppressPackageStartupMessages(require(circlize))

	repTree <- RepeatTree(repRef, method=dendroMethod)
	repDend <- makeBinary(getDendrogram(repTree))
	# repDend <- setMemberAttr(repDend, "label", getMemberAttr(repDend, "label"), unsetInternalNodes=TRUE)

	if (class(compInfo) == "comparisonInfo")   compInfo <- list(compInfo)
	if (is.data.frame(diffScores) || !is.list(diffScores))  diffScores <- list(diffScores)

	nComps <- length(compInfo)
	if (length(diffScores) != nComps) stop(paste0("Invalid argument: diffScores. Must be of length N_comps"))

	if (!is.list(sampleScores))  sampleScores <- list(sampleScores)
	if (length(sampleScores)==1) sampleScores <- rep(sampleScores, nComps)
	if (!is.element(length(sampleScores), c(0, nComps))) stop(paste0("Invalid argument: sampleScores. Must be of length 0, 1 or N_comps"))
	includeSampleScores <- length(sampleScores) > 0

	# convert arguments that can differ for each group to group lists
	if (!is.list(colorGradient.groupScore))  colorGradient.groupScore <- list(colorGradient.groupScore)
	if (length(colorGradient.groupScore)==1) colorGradient.groupScore <- rep(colorGradient.groupScore, nComps)
	if (length(colorGradient.groupScore) != nComps) stop(paste0("Invalid argument: colorGradient.groupScore. Must be of length 1 or N_comps"))

	if (!is.list(zlim.groupScore))  zlim.groupScore <- list(zlim.groupScore)
	if (length(zlim.groupScore)==1) zlim.groupScore <- rep(zlim.groupScore, nComps)
	if (length(zlim.groupScore) != nComps) stop(paste0("Invalid argument: zlim.groupScore. Must be of length 1 or N_comps"))

	if (!is.list(colorGradient.diff))  colorGradient.diff <- list(colorGradient.diff)
	if (length(colorGradient.diff)==1) colorGradient.diff <- rep(colorGradient.diff, nComps)
	if (length(colorGradient.diff) != nComps) stop(paste0("Invalid argument: colorGradient.diff. Must be of length 1 or N_comps"))

	if (!is.list(zlim.diff))  zlim.diff <- list(zlim.diff)
	if (length(zlim.diff)==1) zlim.diff <- rep(zlim.diff, nComps)
	if (length(zlim.diff) != nComps) stop(paste0("Invalid argument: zlim.diff. Must be of length 1 or N_comps"))

	if (!is.list(colorGradient.score))  colorGradient.score <- list(colorGradient.score)
	if (length(colorGradient.score)==1) colorGradient.score <- rep(colorGradient.score, nComps)
	if (length(colorGradient.score) != nComps) stop(paste0("Invalid argument: colorGradient.score. Must be of length 1 or N_comps"))

	if (!is.list(zlim.score))  zlim.score <- list(zlim.score)
	if (length(zlim.score)==1) zlim.score <- rep(zlim.score, nComps)
	if (length(zlim.score) != nComps) stop(paste0("Invalid argument: zlim.score. Must be of length 1 or N_comps"))



	repIds <- getRepeatIds(repRef)
	repIds.unnamed <- repIds
	names(repIds) <- repIds #for the ComplexHeatmap row names
	nReps <- length(repIds)
	if (is.null(leafColors)) leafColors <- rainbow(nReps) # rep("grey", nReps)
	leafColors <- as.matrix(leafColors)
	rownames(leafColors) <- repIds

	# if (getConfigElement("debug")) print(paste0("[DEBUG:]  Constructing repeat tree "))
	treeHm <- Heatmap(leafColors[,1,drop=FALSE], rect_gp=gpar(type = "none"),
		cell_fun = function(j, i, x, y, w, h, fill) {
			grid.circle(x=x, y=y, r=h*60, gp=gpar(col=NA, fill = leafColors[i,1]))
		},
		cluster_rows=repDend, cluster_columns=FALSE,
		show_row_names=FALSE, show_column_names=FALSE, show_heatmap_legend=FALSE,
		name="tree"
		# width=unit(0.1, "npc")
	) + rowAnnotation(labels= anno_text(repIds, which = "row", just=c("left", "center")), width=unit(2, "cm"))#width=unit(0.2, "npc")) #interestingly, for custom dendrograms, repIds must be named
	chm <- treeHm

	for (k in 1:nComps){
		# if (getConfigElement("debug")) print(paste0("[DEBUG:]  adding data for comparison", i))
		X.groupScores <- as.matrix(diffScores[[k]][, c("score.g1", "score.g2")])
		colnames(X.groupScores) <- c(compInfo[[k]]$name.grp1, compInfo[[k]]$name.grp2)
		rownames(X.groupScores) <- repIds.unnamed

		zlimV <- zlim.groupScore[[k]]
		if (is.null(zlimV[[k]])) zlimV <- range(X.groupScores, na.rm=TRUE)
		colR.groupScores <- colorRamp2(seq(zlimV[1], zlimV[2], length.out=length(colorGradient.groupScore[[k]])), colorGradient.groupScore[[k]])

		# if (getConfigElement("debug")) print(paste0("[DEBUG:]    CHK1"))
		X.diff <- as.matrix(diffScores[[k]][, c("diffScore")])
		colnames(X.diff) <- paste0("diff_", compInfo[[k]]$cmpName)
		rownames(X.diff) <- repIds.unnamed

		# if (getConfigElement("debug")) print(paste0("[DEBUG:]    CHK2"))

		zlimV <- zlim.diff[[k]]
		if (is.null(zlimV[[k]])) zlimV <- range(X.diff, na.rm=TRUE)
		colR.diff <- colorRamp2(seq(zlimV[1], zlimV[2], length.out=length(colorGradient.diff[[k]])), colorGradient.diff[[k]])

		colR.scores <- NULL
		Xs <- NULL
		if (includeSampleScores){
			Xs <- sampleScores[[k]][repIds.unnamed, ]
			zlimV <- zlim.score[[k]]
			if (is.null(zlimV[[k]])) zlimV <- range(Xs, na.rm=TRUE)
			colR.scores <- colorRamp2(seq(zlimV[1], zlimV[2], length.out=length(colorGradient.score[[k]])), colorGradient.score[[k]])
		}

		if (includeSampleScores){
			chm <- chm + Heatmap(
				Xs[, compInfo[[k]]$sampleIdx.grp1, drop=FALSE],
				col = colR.scores,
				show_row_names=FALSE, cluster_columns=FALSE,
				name=paste0("score_", compInfo[[k]]$cmpName, "_", compInfo[[k]]$name.grp1)
			)

			chm <- chm + Heatmap(
				Xs[, compInfo[[k]]$sampleIdx.grp2, drop=FALSE],
				col = colR.scores,
				show_row_names=FALSE, cluster_columns=FALSE,
				name=paste0("score_", compInfo[[k]]$cmpName, "_", compInfo[[k]]$name.grp2)
			)
		}

		# chm <- chm + Heatmap(
		# 	X.groupScores[,1,drop=FALSE],
		# 	col = colR.groupScores,
		# 	show_row_names=FALSE,
		# 	name=paste0("gScore_hm_", compInfo[[k]]$cmpName, "_", compInfo[[k]]$name.grp1)
		# )

		chm <- chm + Heatmap(
			X.groupScores,
			col = colR.groupScores,
			show_row_names=FALSE, cluster_columns=FALSE,
			name=paste0("groupScore_", compInfo[[k]]$cmpName)
		)

		# chm <- chm + Heatmap(
		# 	X.diff,
		# 	col = colR.diff,
		# 	show_row_names=FALSE,
		# 	name=paste0("diff_hm_", compInfo[[k]]$cmpName)
		# )

		pValVec <- diffScores[[k]]$diffPval.adj
		pValVec[is.na(pValVec)] <- 1
		names(pValVec) <- repIds.unnamed

		chm <- chm + Heatmap(
			X.diff,
			col=colR.diff,
			cell_fun = function(j, i, x, y, w, h, fill) {
				fillCol <- NA
			    if (pValVec[i] < 0.1) {
					fillCol <- "black"
					gt <- "*"
					if (pValVec[i] < 0.05) gt <- "**"
					if (pValVec[i] < 0.01) gt <- "***"
					grid.text(gt, x, y)
					grid.rect(x, y, w, h, gp = gpar(fill=NA, col=fillCol))
				}
			},
			show_row_names=FALSE, cluster_columns=FALSE,
			name=paste0("diff_hm_", compInfo[[k]]$cmpName)
		)
			

		if (includeSampleScores){
			# chm <- chm + rowAnnotation(
			# 	g1 = row_anno_boxplot(Xs[, compInfo[[k]]$sampleIdx.grp1, drop=FALSE], axis=TRUE, gp=gpar(col="#1b7837")),
			# 	g2 = row_anno_boxplot(Xs[, compInfo[[k]]$sampleIdx.grp2, drop=FALSE], axis=TRUE, gp=gpar(col="#762a83")),
			# 	width=unit(0.2, "npc")
			# )
			
			rg <- range(Xs[, c(compInfo[[k]]$sampleIdx.grp1, compInfo[[k]]$sampleIdx.grp2)], na.rm=TRUE)
			anno_multiple_boxplot <- function(index) {
				pushViewport(viewport(xscale=rg, yscale=c(0.5, nReps+0.5)))
				for(i in index) {
					grid.boxplot(Xs[i, compInfo[[k]]$sampleIdx.grp1], pos=nReps-i+1+0.2, box_width=0.3, gp=gpar(col="#1b7837"), direction="horizontal")
					grid.boxplot(Xs[i, compInfo[[k]]$sampleIdx.grp2], pos=nReps-i+1-0.2, box_width=0.3, gp=gpar(col="#762a83"), direction="horizontal")
				}
				popViewport()
			}
			# chm <- chm + rowAnnotation(boxplot=anno_multiple_boxplot, width=unit(0.2, "npc"), show_annotation_name=FALSE)
			chm <- chm + rowAnnotation(boxplot=anno_multiple_boxplot, width=unit(4, "cm"), show_annotation_name=FALSE)
		}
	}

	# pdftemp(width=20, height=150)
	# 	draw(chm)
	# dev.off()
	return(chm)
}


#' createRepPlot_differential
#'
#' Plots a repeat differential summary tree across all marks for a given comparison in the dataset
#'
#' @param .obj	            \code{\linkS4class{RepeatEpigenomeCollection}} object
#' @param compInfo          comparison info as returned by \code{\link{getComparisonInfo,RepeatEpigenomeCollection-method}}
#' @param plotDir			Output directory where the plots are saved to
#' @param dendroMethod      method for plotting the repeat subfamily dendrogram. See \code{\linkS4class{RepeatTree-class}} class for possible values.
#' @param leafColorMethod   method for coloring the leafs of the dendrogram. Options are \code{"coverage"} for read coverage (default) and \code{"abundance"} for genomic abundance
#' @param normChipMethod    method for normalizing ChIP enrichment scores per mark. See \code{\link{normalizeMatrix}} for possible values. 
#' @param minReads			filtering parameter: specify the minimum number of reads that must match to a given repeat element in order for the repeat to be added to the plot
#' @param minCpGs			filtering parameter: specify the minimum number of CpG that must be contained in a given repeat element in order for the repeat to be added to the plot
#' @param minCpGcov			filtering parameter: specify the minimum number of reads that must cover a given CpG in order to be considered a valid methylation measurement
#' @return nothing of particular interest
#'
#' @details
#' The reference repeats are obtained via \code{RepeatReference()}. Which repeats are used is regulated by the filtering parameters
#' One PDF file is created in the output directory
#'
#' @author Fabian Mueller
#' @noRd
createRepPlot_differential <- function(
		.obj,
		compInfo,
		plotDir,
		dendroMethod="repeatFamily",
		leafColorMethod=getConfigElement("plotRepTree.leafColorMethod"),
		normChipMethod="none",
		minReads=getConfigElement("plotRepTree.meth.minReads"),
		minCpGs=getConfigElement("plotRepTree.meth.minCpGs"),
		minCpGcov=getConfigElement("meth.minCpGcov")){

	rec <- filterRepRefMeth(.obj, minReads=minReads, minCpGs=minCpGs, minCpGcov=minCpGcov)
	rec <- filterRepRefChip(rec, minReads=minReads)
	rec <- filterRepRefAcc(rec, minReads=minReads)
	repRef <- getRepRef(rec)
	if (length(getRepeatIds(repRef)) < 1) logger.error("No repeat sequence retained after filtering")

	markLvls <- getMarks(rec)

	diffMatL <- lapply(markLvls, FUN=function(mn){getRepeatScoresDiff(rec, mn, compInfo)})
	names(diffMatL) <- markLvls
	markIsValid <- sapply(diffMatL, FUN=function(x){!is.null(x)})
	if (sum(markIsValid) < 1) logger.error("No differential matrix could be computed for any mark")
	markLvls <- markLvls[markIsValid]
	diffMatL <- diffMatL[markIsValid]

	scoreMatL <- lapply(markLvls, FUN=function(mn){
		rr <- getRepeatScores(rec, mn, minCpGcov=minCpGcov)
		mt <- inferMarkTypes(mn)
		if (is.element(mt, c("ChIPseq", "Acc"))){
			rr <- normalizeMatrix(rr, method=normChipMethod)
		}
		return(rr)
	})
	names(scoreMatL) <- markLvls

	covgMatL <- lapply(markLvls, FUN=function(mn){
		cm <- getRepeatCovg(rec, mn)
		if (mn == "DNAmeth" && all(is.na(cm))){
			logger.info("No read coverage found. --> using the number of repeat instances as a proxy.")
			cm <- getRepeatCovg(rec, mn, type="numInstances")
		}
		return(cm)
	})
	names(covgMatL) <- markLvls

	leafColScore <- NULL
	if (leafColorMethod=="abundance"){
		repInfo <- getRepeatInfo(repRef)
		leafColScore <- repInfo[,"percCovg"]
	} else {
		# relative coverage
		covgMat.rel <- do.call("cbind", lapply(covgMatL, FUN=function(mm){
			rowMeans(apply(mm, 2, FUN=function(x){x/sum(x)}), na.rm=TRUE)
		})) 
		leafColScore <- rowMeans(covgMat.rel, na.rm=TRUE)
	}
	leafColors <- colorize.value(leafColScore, colscheme=colorpanel(100,"white","black"))
	leafColors[is.na(leafColScore)] <- "#bd0026" #"red"

	colGrads.score <- rep(list(colorpanel(100,"#8C510A","#F5F5F5","#01665E")), length(markLvls))
	colGrads.groupScore <- rep(list(colorpanel(100,"#EDF8B1","#41B6C4","#081D58")), length(markLvls))
	colGrads.diff <- rep(list(colorpanel(100,"#2166ac","#F7F7F7","#b2182b")), length(markLvls))

	zlims.score <- rep(list(c(-3, 3)), length(markLvls)) #chip-seq abs(log2FC) limit to range [-3,3]
	zlims.groupScore <- rep(list(NULL), length(markLvls))
	zlims.diff <- rep(list(c(-4, 4)), length(markLvls))

	indMeth <- match("DNAmeth", markLvls)
	if (!is.na(indMeth)){
		colGrads.score[[indMeth]] <- colorpanel(100,"#EDF8B1","#41B6C4","#081D58")
		zlims.score[[indMeth]] <- c(0,1)
		zlims.diff[[indMeth]] <- c(-1,1)
	}
	indAcc <- match("Acc", inferMarkTypes(markLvls))
	if (!is.na(indAcc) && length(indAcc)){
		for (i in indAcc){
			zlims.score[[i]] <- c(-5,5)
		}
	}

	chm <- repPlot_differential(
		repRef,
		compInfo,
		diffMatL,
		sampleScores=scoreMatL,
		colorGradient.groupScore=colGrads.groupScore,
		zlim.groupScore=zlims.groupScore,
		colorGradient.diff=colGrads.diff,
		zlim.diff=zlims.diff,
		colorGradient.score=colGrads.score,
		zlim.score=zlims.score,
		leafColors=leafColors,
		dendroMethod=dendroMethod
	)

	fn <- file.path(plotDir, paste0("repeatTree_differential_", compInfo$cmpName, ".pdf"))
	pdf(fn, width=20, height=150)
		draw(chm)
	dev.off()

	res <- list(
		status="success",
		callParams=list(
			plotDir=plotDir,
			minReads=minReads,
			minCpGs=minCpGs,
			minCpGcov=minCpGcov
		)
	)
}
