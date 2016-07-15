################################################################################
# methods for fasta-based annotation of repeats from RepBaseUpdate
################################################################################

#' getDNAStringForReference
#'
#' given a reference FASTA file, retrieve sequence information for each RE
#'
#' @param reference Filename of the reference of REs (FASTA file)
#' @return \code{DNAStringSet} object containing sequence information for each RE
#'
#' @author Fabian Mueller
#' @noRd
getDNAStringForReference <- function(reference){
	refSeqs <- readBStringSet(reference)
	#replace letters "x"
	for (i in 1:length(refSeqs)){
		refSeqs[[i]][start(matchPattern("x",refSeqs[[i]]))] <- "n"
	}
	refSeqs <- DNAStringSet(refSeqs)
	return(refSeqs)
}
#' getReferenceInfo
#'
#' given a reference FASTA file, retrieve sequence annotation for each RE
#'
#' @param reference Filename of the reference of REs (FASTA file)
#' @return \code{data.frame} containing annotation for each RE
#' 
#' @details the RE identifiers in the FASTA file are assumed to have the following format
#' (tab separated):
#' repeat id, repeat family, species
#'
#' @author Fabian Mueller
#' @noRd
getReferenceInfo <- function(reference){
	refSeqs <- readBStringSet(reference)
	seqLens <- vapply(refSeqs,length,integer(1))
	res <- data.frame(t(data.frame(strsplit(names(refSeqs),"\t"))),seqLength=seqLens)
	colnames(res)[1:4] <- c("id","family","species","seqLength")
	res[,"id"] <- as.character(res[,"id"])
	rownames(res) <- res[,1]
	return(res)
}

################################################################################
# Methods for manual data curation
################################################################################

#' getTaxOrder
#'
#' gets a taxonomical ordering for a given taxonomical term
#' from the hand-curated annotation
#' 
#' @param taxVec   a character or factor vector of taxonomical terms which
#'                 should be put in order
#' @param tax      target taxonomy/species for the simplification scheme
#'
#' @return a factor vector with levels ordered according to the hand-curated annotation
#'
#' @author Fabian Mueller
#' @noRd
getTaxOrder <- function(taxVec, tax="human"){
	taxVec <- as.character(taxVec)
	cureScheme <- fromJSON(system.file(file.path("extdata", "curation", "taxonomyAnnotation.json"), package="epiRepeatR"))
	if (!is.element(tax, names(cureScheme))){
		logger.error(c("taxonomy term not found in curation scheme:",tax))
	}
	cur <- cureScheme[[tax]][["taxOrder"]]
	res <- factor(taxVec, levels=cur)
	return(res)
}

#' simplifyRepeatFamilies
#'
#' simplifies the annotated repeat families to get rid of annotation artifacts
#' 
#' @param famVec   a character or factor vector of repeat families
#' @param tax      taxonomy/species for the simplification scheme
#'
#' @return a factor vector containing the simplified repeat families
#'
#' @author Fabian Mueller
#' @noRd
simplifyRepeatFamilies <- function(famVec, tax="human"){
	famVec <- as.character(famVec)
	cureScheme <- fromJSON(system.file(file.path("extdata", "curation", "repeatFamilies.json"), package="epiRepeatR"))
	if (!is.element(tax, names(cureScheme))){
		logger.error(c("taxonomy term not found in curation scheme:",tax))
	}
	res <- famVec
	cur <- cureScheme[[tax]][["familyRegExp"]]
	curFams <- names(cur)
	for (fn in curFams){
		for (i in 1:length(cur[[fn]])){
			curPat <- paste0("^",cur[[fn]][i],"$")
			res[grepl(curPat, famVec)] <- fn
		}
	}
	res <- factor(res)
	return(res)
}

#' getCuratedRepeatFamilyTree
#'
#' retrieves a hand-curated family tree of repeat families
#'
#' @return a dendrogram containing the curated repeat family tree
#'
#' @author Fabian Mueller
#' @noRd
getCuratedRepeatFamilyTree <- function(){
	treeList <- fromJSON(system.file(file.path("extdata", "curation", "repeatTree.json"), package="epiRepeatR"))
	if (length(treeList) > 1){
		logger.error("Not a tree (but a forest?)")
	}
	treeName <- "Repetitive element"
	res <- list2dend(treeList[[treeName]], treeName)
	res <- adjustAttr.midpoint(res)
	return(res)
}
# rtc <- getCuratedRepeatFamilyTree()
# plot(rtc)
# require(diagram)
# openplotmat()
# plotDend.rec(rtc, xmin=0, xmax=1, ymin=0, ymax=1, rev=TRUE)

#' addLeafs
#'
#' create a tree that contains all the child nodes of the original tree
#' plus additional leaf nodes for a vector of new node labels
#'
#' @param dend    dendrogram to modify
#' @param newLeafLabels labels for the newly created leaf nodes.
#'                if empty, the unmodified dendrogram is returned
#' 
#' @return the modified dendrogram
#' 
#' @author Fabian Mueller
#' @noRd
addLeafs <- function(dend, newLeafLabels){
	res <- dend
	if (length(newLeafLabels) > 0){
		dends2merge <- lapply(newLeafLabels, createStump)
		# the dendrogram is not a leaf anymore
		if (!isLeaf(dend)){
			subTrees <- lapply(dend, identity)
			dends2merge <- c(subTrees, dends2merge)
		}
		if (length(dends2merge) > 1){
			# print(sapply(dends2merge,class))
			res <- do.call("merge", unname(dends2merge))
		} else if (length(dends2merge) == 1){
			res <- createLinearTree(dends2merge[[1]], attr(dend, "label"), h=attr(dend, "height"))
		} else {
			logger.error("Invalid structure: cannot combine less than one node into a tree")
		}
		nMems <- attr(res, "members")
		attr(res, "label")    <- attr(dend, "label")
		attr(res, "height")   <- attr(dend, "height")
		if (isLeaf(dend)) attr(res, "leaf") <- NULL
		attr(res, "midpoint") <- (nMems-1)/2
	}
	return(res)
}
#' assembleRepeatsInCuratedFamilyTree
#'
#' place the given repeats into the curated family tree
#' given their annotated family
#'
#' @param ids     character vector of repeat ids
#' @param fams    character vector of repeat families.
#'                must be matching and of same length as \code{ids}
#' @param famDend dendrogram of repeat families
#' 
#' @return a dendrogram containing the curated repeat family tree with the
#'         repeat ids as leafs
#' 
#' @details ! Recursive !
#' 
#' @author Fabian Mueller
#' @noRd
assembleRepeatsInCuratedFamilyTree <- function(ids, fams, famDend=getCuratedRepeatFamilyTree()){
	# print(attr(famDend,"label"))
	matchingFams <- fams==attr(famDend,"label")
	curIds <- c()
	if (any(matchingFams)){
		curIds <- ids[matchingFams]
	}
	res <- famDend
	if (!isLeaf(famDend)){
		dendList <- lapply(famDend, FUN=function(x){
			assembleRepeatsInCuratedFamilyTree(ids, fams, famDend=x)
		})
		if (length(dendList) > 1){
			res <- do.call("merge", unname(dendList))
		} else if (length(dendList) == 1){
			res <- createLinearTree(dendList[[1]], attr(famDend,"label"), attr(famDend,"height"))
		} else {
			logger.error("Invalid structure: no subTrees found and not a leaf")
		}
	}
	res <- addLeafs(res, curIds)
	attr(res, "height") <- attr(famDend, "height") + 1L
	return(res)
}
# rt <- assembleRepeatsInCuratedFamilyTree(repFeats[["id"]], simplifyRepeatFamilies(repFeats[["family"]], tax="human"))
# rt <- adjustAttr.midpoint(rt)
# pdf("~/tmp/repTreeCure_dendPlot.pdf", width=100,height=10)
# 	plot(rt)
# dev.off()
# require(diagram)
# pdf("~/tmp/repTreeCure.pdf", width=10,height=100)
# 	openplotmat()
# 	plotDend.rec(rt, xmin=0, xmax=1, ymin=0, ymax=1, rev=TRUE)
# dev.off()
