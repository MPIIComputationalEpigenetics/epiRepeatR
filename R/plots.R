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
	image(
		gridlinePos.x, gridlinePos.y,
		t(X),
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
		addBoxplots=FALSE){

	require(diagram)
	repTree <- RepeatTree(repRef, method="repeatFamily")
	# repTree <- RepeatTree(repRef, method="hierClust")
	orderedLabels <- getDendrogramMembers(repTree, rev=FALSE)
	orderedLabels4plot <- getDendrogramMembers(repTree, rev=TRUE)
	scores.ordered4plot <- scores[orderedLabels4plot,]
	
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
		text(heatResMean$gridMidPoints.x, coord.heat.grp.sum["ymax",i], summaryColNames, srt=45,adj=c(0,0), cex=txt.cex.base, col=groupColors[i])

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

filterRepRef_meth <- function(repRef, methCallResList, minReads=100, minCpGs=2){
	repRefNames <- getRepeatIds(repRef)

	# rrKmerCounts <- getKmerCounts(repRef)
	numCG <- do.call("cbind",lapply(methCallResList, FUN=function(x){
		sapply(x[repRefNames], FUN=function(r){
			if (is.null(r)) return(0) else return(length(r$methCalls$cPos))
		})
	}))
	numReads <- do.call("cbind",lapply(methCallResList, FUN=function(x){
		sapply(x[repRefNames], FUN=function(r){
			if (is.null(r)) return(0) else return(r$readStats["numReads"])
		})
	}))

	survive <- numCG >= minCpGs & numReads >= minReads
	survive <- apply(survive,1,all)
	return(filterRepeats_wl(repRef, repRefNames[survive]))
}
filterRepRef_chip <- function(repRef, quantResList, minReads=100){
	repRefNames <- getRepeatIds(repRef)

	numReads.chip <- do.call("cbind",lapply(quantResList, FUN=function(x){
		sapply(x[repRefNames], FUN=function(r){
			if (is.null(r)) return(0) else return(r$readStats["numReads_chip"])
		})
	}))
	numReads.input <- do.call("cbind",lapply(quantResList, FUN=function(x){
		sapply(x[repRefNames], FUN=function(r){
			if (is.null(r)) return(0) else return(r$readStats["numReads_input"])
		})
	}))

	survive <- numReads.chip >= minReads & numReads.input >= minReads
	survive <- apply(survive,1,all)
	return(filterRepeats_wl(repRef, repRefNames[survive]))
}

#' createRepPlot_groupSummaryTrees_meth
#'
#' Plots a repeat sequence group summary tree for methylation values for each grouping information entry
#'
#' @param repMethCallFns	vector of filenames which point to RData objects (\code{rds} files) containing methylation information. They should contain object in the format as ouput by \code{\link{getMethylationCalls,RepeatAlignment-method}} 
#' @param sampleNames		vector of sample names to be used. Must be in the same order as \code{repMethCallFns}
#' @param sampleGroups 		list of assignments of samples/experiments to groups. Obtained by the \code{\link{getSampleGroups}} function
#' @param plotDir			Output directory where the plots are saved to
#' @param minReads			filtering parameter: specify the minimum number of reads that must match to a given repeat element in order for the repeat to be added to the plot
#' @param minCpGs			filtering parameter: specify the minimum number of CpG that must be contained in a given repeat element in order for the repeat to be added to the plot
#' @return nothing of particular interest
#'
#' @details
#' The reference repeats are obtained via \code{RepeatReference()}. Which repeats are used is regulated by the filtering parameters
#' One PDF file for each comparison information is created in the output directory
#'
#' @author Fabian Mueller
#' @noRd
createRepPlot_groupSummaryTrees_meth <- function(repMethCallFns, sampleNames, sampleGroups, plotDir, minReads=100, minCpGs=2){

						# ft <- getFileTable(am)
						# ft.inds <- ft[,"mark"]=="DNAmeth" & ft[,"fileType"]=="rds" & ft[,"analysisStep"]=="methCalling"
						# repMethCallFns <- prependAnalysisDirForFilename(am, ft[ft.inds,"fileName"], outDir)
						# sampleNames <- ft[ft.inds, "sampleName"]
						# samAnnot <- getSampleAnnot(am)
						# sampleGroups <- getSampleGroups(samAnnot[rownames(samAnnot) %in% sampleNames,,drop=FALSE])
						# plotDir <- "/DEEP_fhgfs/projects/fmueller/tmp/repeatPlots" #file.path(outDir)
						# minReads <- 100
						# minCpGs <- 2


	if (length(repMethCallFns)!=length(sampleNames)){
		stop("incompatible repMethCallFns and sampleNames parameters")
	}

	repRef <- RepeatReference()

	methCallResList <- lapply(repMethCallFns, FUN=function(fn){
		res <- readRDS(fn)
		return(res)
	})
	names(methCallResList) <- sampleNames

	repRef <- filterRepRef_meth(repRef, methCallResList, minReads=minReads, minCpGs=minCpGs)

	repRefNames <- getRepeatIds(repRef)

	methScores <- do.call("cbind",lapply(methCallResList,FUN=function(mc){
		sapply(mc[repRefNames],FUN=function(r){
			if (is.null(r)) {
				return(NA)
			} else {
				methLvl <- mean(r$methCalls[,"numM"]/r$methCalls[,"numT"],na.rm=TRUE)
				return(methLvl)
			}
		})
	}))
	covgMat <- do.call("cbind",lapply(methCallResList,FUN=function(mc){
		sapply(mc[repRefNames],FUN=function(r){
			if (is.null(r)) {
				return(NA)
			} else {
				r$readStats["numReads"]
			}
		})
	}))
	rownames(covgMat) <- repRefNames

	#create a plot for each grouping info
	for (i in 1:length(sampleGroups)){
		ggn <- names(sampleGroups)[i]
		ggs <- sampleGroups[[i]]

		covgMat.sub <- covgMat[,unlist(ggs)]
		covgMat.sub.rel <-  apply(covgMat.sub, 2, FUN=function(x){x/sum(x)}) # relative coverage
		covg.score <- rowMeans(covgMat.sub.rel, na.rm=TRUE)
		#order the color values
		leafColors <- colorize.value(covg.score, colscheme=colorpanel(100,"white","black"))
		fn <- file.path(plotDir, paste0("repeatTree_groupSummary_",ggn,".pdf"))
		# load_all(file.path(pkg.dir,"epiRepeatR"))
		pdf(fn, width=20, height=100)
			repPlot_groupSummary(repRef, methScores, ggs, zlim=c(0,1), leafColors=leafColors, addBoxplots=TRUE)
		dev.off()
	}
	res <- list(
		status="success",
		callParams=list(
			repMethCallFns=repMethCallFns,
			sampleNames=sampleNames,
			sampleGroups=sampleGroups,
			plotDir=plotDir,			
			minReads=minReads,
			minCpGs=minCpGs
		)
	)
}


#' createRepPlot_markTree
#'
#' Plots a repeat sequence summary tree for all marks in the dataset
#'
#' @param repMethCallFns	vector of filenames which point to RData objects (\code{rds} files) containing methylation information. They should contain object in the format as ouput by \code{\link{getMethylationCalls,RepeatAlignment-method}} 
#' @param sampleNames		vector of sample names to be used. Must be in the same order as \code{repMethCallFns}
#' @param markNames 		vector of names of the mark assayed
#' @param plotDir			Output directory where the plots are saved to
#' @param minReads			filtering parameter: specify the minimum number of reads that must match to a given repeat element in order for the repeat to be added to the plot
#' @param minCpGs			filtering parameter: specify the minimum number of CpG that must be contained in a given repeat element in order for the repeat to be added to the plot
#' @return nothing of particular interest
#'
#' @details
#' The reference repeats are obtained via \code{RepeatReference()}. Which repeats are used is regulated by the filtering parameters
#' One PDF file for each comparison information is created in the output directory
#'
#' @author Fabian Mueller
#' @noRd
createRepPlot_markTree <- function(quantFns, sampleNames, markNames, plotDir, minReads=100, minCpGs=2){

	if (length(quantFns)!=length(sampleNames) || length(quantFns)!=length(markNames)) {
		stop("incompatible quantFns, sampleNames and markNames parameters")
	}
	sampleNames.full <- paste(sampleNames, markNames, sep="_")

	repRef <- RepeatReference()
	nFiles <- length(quantFns)
	quantResList <- lapply(quantFns, FUN=function(fn){
		res <- readRDS(fn)
		return(res)
	})
	names(quantResList) <- sampleNames.full

	isMeth <- markNames=="DNAmeth"
	isChip <- !isMeth

	repRef <- filterRepRef_meth(repRef, quantResList[isMeth], minReads=minReads, minCpGs=minCpGs)
	repRef <- filterRepRef_chip(repRef, quantResList[isChip], minReads=minReads)

	repRefNames <- getRepeatIds(repRef)

	scoreMat <- do.call("cbind",lapply(1:nFiles, FUN=function(i){
		sapply(quantResList[[i]][repRefNames],FUN=function(r){
			if (is.null(r)) {
				return(NA)
			} else {
				if (isMeth[i]){
					sc <- mean(r$methCalls[,"numM"]/r$methCalls[,"numT"],na.rm=TRUE)
				} else if (isChip[i]){
					sc <- r$log2fc
				} else {
					sc <- NA
				}
				return(sc)
			}
		})
	}))
	colnames(scoreMat) <- sampleNames.full

	#color leafs by mean relative coverage
	covgMat <- do.call("cbind",lapply(1:nFiles, FUN=function(i){
		sapply(quantResList[[i]][repRefNames],FUN=function(r){
			if (is.null(r)) {
				return(NA)
			} else {
				if (isMeth[i]){
					res <- r$readStats["numReads"]
				} else if (isChip[i]){
					res <- r$readStats["numReads_chip"]
				} else {
					res <- NA
				}
				return(res)
			}
		})
	}))
	colnames(covgMat) <- sampleNames.full
	covgMat.rel <- apply(covgMat, 2, FUN=function(x){x/sum(x)}) # relative coverage
	covg.score <- rowMeans(covgMat.rel, na.rm=TRUE)
	leafColors <- colorize.value(covg.score, colscheme=colorpanel(100,"white","black"))

	markLvls <- sort(unique(markNames))
	ggs <- lapply(markLvls, FUN=function(mn){
		sampleNames.full[markNames==mn]
	})
	names(ggs) <- markLvls

	colGrads <- rep(list(colorpanel(100,"#8C510A","#F5F5F5","#01665E")), length(markLvls))
	zlims <- rep(list(c(-3, 3)), length(markLvls)) #chip-seq abs(log2FC) limit to range [-3,3]
	indMeth <- match("DNAmeth", markLvls)
	if (!is.na(indMeth)){
		colGrads[[indMeth]] <- colorpanel(100,"#EDF8B1","#41B6C4","#081D58")
		zlims[[indMeth]] <- c(0,1)
	}

	fn <- file.path(plotDir, paste0("repeatTree_markSummary",".pdf"))
	pdf(fn, width=20, height=100)
		repPlot_groupSummary(repRef, scoreMat, ggs, colorGradient=colGrads, zlim=zlims, leafColors=leafColors, addBoxplots=FALSE)
	dev.off()

	res <- list(
		status="success",
		callParams=list(
			quantFns=quantFns,
			sampleNames=sampleNames,
			markNames=markNames,
			plotDir=plotDir,			
			minReads=minReads,
			minCpGs=minCpGs
		)
	)
}

