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
#' recursive function transforming a list to a dendrogram
#'
#' @return a dendrogram
#' 
#' @details Leaf nodes are represented by \code{list()}
#'
#' @author Fabian Mueller
#' @noRd
list2dend <- function(x, lab="[root]"){
	res <- NULL
	if (length(x)<1){
		res <- lab
		attr(res,"members") <- 1L
		attr(res,"label") <- lab
		attr(res,"height") <- 0L
		attr(res,"leaf") <- TRUE
	} else if (length(x)==1){
		res <- list(list2dend(x[[1]], names(x)[1]))
		attr(res,"members")  <- attr(res[[1]],"members")
		attr(res,"height")   <- attr(res[[1]],"height") + 1L
		# attr(res,"midpoint") <- 1
		attr(res,"label")   <- lab
	} else {
		subTrees <- lapply(names(x), FUN=function(ll){list2dend(x[[ll]], ll)})
		stMembers <- vapply(subTrees, FUN=function(y){attr(y,"members")}, integer(1))
		stHeights <- vapply(subTrees, FUN=function(y){attr(y,"height")}, integer(1))

		res <- do.call("merge", unname(subTrees))

		attr(res,"members")  <- sum(stMembers)
		attr(res,"label")    <- lab
		attr(res,"height")   <- max(stHeights) + 1L
		attr(res,"midpoint") <- (attr(res,"members")-1)/2
	}
	class(res) <- "dendrogram"
	return(res)
}
# dd <- list2dend(treeList)
# plot(dd)

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
	return(list2dend(treeList, names(treeList)[1]))
}

# #' assembleRepeatsInCuratedFamilyTree
# #'
# #' place the given repeats into the curated family tree
# #' given their annotated family
# #'
# #' @param ids     character vector of repeat ids
# #' @param fams    character vector of repeat families.
# #'                must be matching and of same length as \code{ids}
# #' @param famDend dendrogram of repeat families
# #' 
# #' @return a dendrogram containing the curated repeat family tree with the
# #'         repeat ids as leafs
# #' 
# #' @details ! Recursive !
# #' 
# #' @author Fabian Mueller
# #' @noRd
# assembleRepeatsInCuratedFamilyTree <- function(ids, fams, famDend=getCuratedRepeatFamilyTree()){
# 	res <- NULL
# 	if (isLeaf(famDend)){
# 		a <- NULL
# 	} else {
# 		a <- NULL
# 	}
# 	return(res)
# }
# assembleRepeatsInCuratedFamilyTree(repFeats[["id"]], simplifyRepeatFamilies(repFeats[["family"]], tax="human"))
