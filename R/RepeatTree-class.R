setOldClass(c("dendrogram"), prototype=structure(list(), class="dendrogram"))
setClassUnion("ListOrNULL", c("list", "NULL"))

#' RepeatTree Class
#'
#' A class for managing a collection of reference repeats in a dendrogram structure
#'
#' @section Slots:
#' \describe{
#'   \item{\code{tree}}{the dendrogram organizing the reference repeats}
#'   \item{\code{repeatRef}}{Reference repeats as Object of type \code{\linkS4class{RepeatReference}}} 
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{\code{\link{getDendrogram,RepeatTree-method}}}{retrieve the repeat dendrogram}
#' }
#'
#' @noRd
#' @name RepeatTree-class
#' @rdname RepeatTree-class
#' @author Fabian Mueller
setClass("RepeatTree",
	slots = list(
		tree="dendrogram",
		repeatRef="RepeatReference"
	),
	package = "epiRepeatR"
)
setMethod("initialize", "RepeatTree",
	function(.Object,
		repRef,
		method="repeatFamily"
	){
		repTree <- NULL
		kmerCounts <- getKmerCounts(repRef)
		dd <- dist(t(kmerCounts), method="euclidean")
		repInfo <- getRepeatInfo(repRef)
		if (method=="repeatFamily"){
			fams <- sort(as.character(unique(repInfo$family)))
			repTree <- lapply(fams,FUN=function(rf){
				iis <- which(repInfo$family==rf)
				res <- lapply(iis,FUN=function(i){
					ress <- i
					attr(ress,"members") <- 1L
					attr(ress,"label") <- repInfo$id[i]
					attr(ress,"height") <- 0
					attr(ress,"leaf") <- TRUE
					class(ress) <- "dendrogram"
					return(ress)
				})
				attr(res,"members") <- length(iis)
				attr(res,"label") <- rf
				attr(res,"height") <- 0.5
				attr(res,"midpoint") <- (length(iis)-1)/2
				class(res) <- "dendrogram"
				return(res)
			})
			attr(repTree,"members") <- nrow(repInfo)
			attr(repTree,"height") <- 1
			attr(repTree,"midpoint") <- (nrow(repInfo)-1)/2
			class(repTree) <- "dendrogram"
		} else if (method=="hierClust"){
			clustHier <- hclust(dd, method="complete")
			repTree <- as.dendrogram(clustHier)
			# repTree <- unclass(repTree)

		# } else if (method=="nj"){
		# 	require(ape)
		# 	phyTree <- nj(dd)
		# 	# pp <- plot(phyTree, type="p", font=1, cex=0.1, use.edge.length=FALSE, show.node.label=TRUE)
		# 	# phyTree$edge.length <- rep(0,length(phyTree$edge.length)) #make ultrametric
		# 	# phyTree$root.edge <- 0 #make the tree rooted
		# 	clustHier <- as.hclust(phyTree) #doesn't work: tree is not ultrametric, rooted or binary
		# 	# TODO: implement further...
		} else if (method=="annotClust"){
			if (is.null(repRef@repeatInfoList)){
				logger.info("Adding repeat annotation from EMBL file...")
				repRef <- addRepeatInfoFromEmbl(repRef)
			}
			ril <- repRef@repeatInfoList
			repTerms <- unique(unlist(lapply(ril, FUN=function(x){x$repeatTerms})))
			repTermMat <- t(sapply(ril, FUN=function(x){
				repTerms %in% x$repeatTerms
			}))
			colnames(repTermMat) <- repTerms
			rownames(repTermMat) <- names(ril)

			termDist <- dist(repTermMat, method="manhattan")
			clustHier <- hclust(termDist, method="single")
			repTree <- as.dendrogram(clustHier)
		} else {
			stop("Unknown method")
		}
		.Object@tree=repTree
		.Object@repeatRef=repRef
		.Object
	}
)
#' @param repRef    Object of type \code{\linkS4class{RepeatReference}}
#' @param method    Method for constructing the dendrogram. Valid methods include 
#'                  grouping by repeat subfamily ("repeatFamily"),
#'                  grouping by hierarchical clustering based on k-mer counts in the repeat sequence (Euclidean distance, complete linkage) ("hierClust"),
#'                  grouping by hierarchical clustering based on occurrences of terms in the annotation fields of a repeat ("annotClust"; requires that the repeat references has been annotated from the EMBL format.)
#' @name RepeatTree
#' @rdname RepeatTree-class
#' @aliases intialize,RepeatTree-method
#' @noRd
RepeatTree <- function(repRef,method="repeatFamily"){
	obj <- new("RepeatTree",
		repRef,
		method=method
	)
	return(obj)
}

if (!isGeneric("getDendrogram")) setGeneric("getDendrogram", function(.Object) standardGeneric("getDendrogram"))
#' getDendrogram-methods
#'
#' Return the dendrogram of the object
#'
#' @param .Object \code{\linkS4class{RepeatTree}} object
#' @return Dendrogram object containing reference repeats
#'
#' @rdname getDendrogram-RepeatTree-method
#' @docType methods
#' @aliases getDendrogram
#' @aliases getDendrogram,RepeatTree-method
#' @author Fabian Mueller
#' @noRd
setMethod("getDendrogram", signature(.Object="RepeatTree"),
	function(.Object){
		return(.Object@tree)
	}
)

if (!isGeneric("getDendrogramMembers")) setGeneric("getDendrogramMembers", function(.Object, ...) standardGeneric("getDendrogramMembers"))
#' getDendrogramMembers-methods
#'
#' Return the repeat members stored in the dendrogram
#'
#' @param .Object \code{\linkS4class{RepeatTree}} object
#' @param rev     should the names be returned in reverse order
#' @return character vector of reference repeat identifiers
#'
#' @rdname getDendrogramMembers-RepeatTree-method
#' @docType methods
#' @aliases getDendrogramMembers
#' @aliases getDendrogramMembers,RepeatTree-method
#' @author Fabian Mueller
#' @noRd
setMethod("getDendrogramMembers", signature(.Object="RepeatTree"),
	function(.Object, rev=FALSE){
		res <- getMemberAttr(.Object@tree, "label")
		if (rev) res <- rev(res)
		return(res)
	}
)

if (!isGeneric("setDendrogramMemberAttr")) setGeneric("setDendrogramMemberAttr", function(.Object, ...) standardGeneric("setDendrogramMemberAttr"))
#' setDendrogramMemberAttr-methods
#'
#' Sets attributes for the members/leafs of the dendrogram
#'
#' @param .Object  \code{\linkS4class{RepeatTree}} object
#' @param attrName name of the attribute to be set
#' @param attrVals values of the attribute to be set for each member
#' @return the modified object with the modified dendrogram containing the modified attributes
#'
#' @rdname setDendrogramMemberAttr-RepeatTree-method
#' @docType methods
#' @aliases setDendrogramMemberAttr
#' @aliases setDendrogramMemberAttr,RepeatTree-method
#' @author Fabian Mueller
#' @noRd
setMethod("setDendrogramMemberAttr", signature(.Object="RepeatTree"),
	function(.Object, attrName, attrVals){
		.Object@tree <- setMemberAttr(.Object@tree, attrName, attrVals)
		return(.Object)
	}
)

# Helper function determining if the dendrogram is a leaf 
isLeaf <- function(dend){
	if (!is.null(attr(dend,"leaf"))) return(attr(dend,"leaf"))
	return(FALSE)
}

#' plotDend.rec
#'
#' Plots a a dendrogram in the specified set of coordinates. Recursive function!
#'
#' @param dend      dendrogram object
#' @param xmin      lower x position where the dendrogram should be plotted
#' @param xmax      upper x position where the dendrogram should be plotted
#' @param ymin      lower y position where the dendrogram should be plotted
#' @param ymax      upper y position where the dendrogram should be plotted
#' @param depth     the current depths in the dendrogram. Important for recursion
#' @param yroot     the y coordinate of the root of the dendrogram
#' @param rev       should the order of the members of the dendrogram be reversed in plotting (top to bottom instead of bottom to top drawing)
#' @param cex       text size of the labels
#' @return nothing of particular interest
#'
#' @author Fabian Mueller
#' @noRd
plotDend.rec <- function(dend, xmin=0, xmax=1, ymin=0, ymax=1, depth=0, yroot=NULL, rev=TRUE, cex=0.5){
	if (isLeaf(dend)){
		pbgcol <- attr(dend,"leafBgColor")
		if (is.null(pbgcol)){
			# pbgcol <- "white"
			pbgcol <- "black"
		}
		leaf.x <- (xmax-xmin)/2+xmin
		leaf.y <- (ymax-ymin)/2+ymin
		points(leaf.x, leaf.y, pch=21, col="black", bg=pbgcol)
		# textempty(cbind(leaf.x, leaf.y), lab=attr(dend,"label"), cex=cex, adj=c(0,0.5), box.col=pbgcol, col=getFgColorForBg(pbgcol))
		textempty(cbind(leaf.x, leaf.y), lab=attr(dend,"label"), cex=cex, adj=c(0,0.5))
		a <- NULL
	} else {
		# print(dend)
		# print(paste0("xmin: ",xmin))
		# print(paste0("xmax: ",xmax))
		# print(paste0("ymin: ",ymin))
		# print(paste0("ymax: ",ymax))
		# print("########################################################################")
		rect(xmin,ymin,xmax,ymax, border="pink")
		
		Nmem <- attr(dend,"members")
		Nbranches <- length(dend)

		root.x <- xmin # derived from attr(dend,"height")
		if (is.null(yroot)){
			root.y <- attr(dend,"midpoint")
			if (rev) {
				root.y <- (Nmem-1) - root.y
			}
			root.y <- root.y / (Nmem-1) * (ymax-ymin) + ymin
		} else {
			root.y <- yroot
		}

		subDend.n <- sapply(dend,FUN=function(x){attr(x,"members")})
		if (rev) subDend.n <- rev(subDend.n)

		height.prop <- (attr(dend,"height") - sapply(dend,FUN=function(x){attr(x,"height")}))/attr(dend,"height")
		subDend.x <- (xmax-xmin) * height.prop + xmin
		if (rev) subDend.x <- rev(subDend.x)

		subDend.ymin <- ymin
		subDend.ymax <- ymax
		subDend.y <- root.y
		#handle the case where the dendrogram just has 1 subbranch
		if (length(dend)>1){
			subDend.n.cum <- cumsum(subDend.n)
			subDend.n.cum.off <- c(0, subDend.n.cum[1:(length(subDend.n.cum)-1)])
			subDend.ymin <- (subDend.n.cum.off / Nmem) * (ymax-ymin) + ymin
			subDend.ymax <- (subDend.n.cum / Nmem) * (ymax-ymin) + ymin
			mids.offset <- sapply(dend,FUN=function(x){
				if(isLeaf(x) || length(x)<2){
					return(0.5)
				} else {
					mp <- attr(x,"midpoint")
					if (rev) mp <-  (attr(x,"members") - 1) - mp
					return(mp / (attr(x,"members") - 1))
				}
			})
			if (rev) mids.offset <- rev(mids.offset)
			subDend.y <- mids.offset * (subDend.ymax - subDend.ymin) + subDend.ymin
		}


		# splitarrow(from=cbind(root.x,root.y),to=cbind(subDend.x,subDend.y),lcol="red")
		# treearrow(from=cbind(root.x,root.y),to=cbind(subDend.x,subDend.y),path="V",arr.side=0,lcol="red")
		treefork(from=cbind(root.x,root.y),to=cbind(subDend.x,subDend.y),path="H")

		for (i in 1:length(dend)){
			subi <- i
			if (rev) subi <- length(dend) - i + 1
			plotDend.rec(dend[[subi]], xmin=subDend.x[i], xmax=xmax, ymin=subDend.ymin[i], ymax=subDend.ymax[i], depth=depth+1, yroot=subDend.y[i], rev=rev)
		}
	}
}

if (!isGeneric("addToPlot")) setGeneric("addToPlot", function(.Object, ...) standardGeneric("addToPlot"))
#' addToPlot-methods
#'
#' Add a repeat tree dendrogram to a plot
#'
#' @param .Object   \code{\linkS4class{RepeatTree}} object
#' @param xmin      lower x position where the dendrogram should be plotted
#' @param xmax      upper x position where the dendrogram should be plotted
#' @param ymin      lower y position where the dendrogram should be plotted
#' @param ymax      upper y position where the dendrogram should be plotted
#' @param rev       should the order of the members of the dendrogram be reversed in plotting (top to bottom instead of bottom to top drawing)
#' @return nothing of particular interest
#'
#' @rdname addToPlot-RepeatTree-method
#' @docType methods
#' @aliases addToPlot
#' @aliases addToPlot,RepeatTree-method
#' @author Fabian Mueller
#' @noRd
setMethod("addToPlot", signature(.Object="RepeatTree"),
	function(.Object, xmin=0, xmax=1, ymin=0, ymax=1, rev=TRUE){
		require(diagram)
		dend <- .Object@tree
		# openplotmat()
		# xmin=0; xmax=1; ymin=0; ymax=1
		plotDend.rec(dend, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, rev=rev)
	}
)

