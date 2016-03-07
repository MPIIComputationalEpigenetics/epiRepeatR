#' dendrogramMergeSimple
#'
#'  naively merge a list of dendrograms into a new dendrogram with each branch belonging to an element in the list by just setting the class and a couple of attributes
#'
#' @param dendroList	list of dendrograms
#' @param height		height (will be set as attribute) on which to merge the dendrograms
#' @param midpoint		midpoint attribute that will be set in the resulting dendrogram
#' @return a dendrogram object with the dendrograms from the input list as branches
#'
#' @author Fabian Mueller
#' @noRd
dendrogramMergeSimple <- function(dendroList, height=NA, midpoint=NULL){
	res <- dendroList
	subTreeMems <- sapply(dend,FUN=function(x){attr(x,"members")})
	attr(res,"members") <- sum(subTreeMems)
	attr(res,"height") <- height
	if (is.null(midpoint)){
		attr(res,"midpoint") <- (sum(subTreeMems)-1)/2
	} else {
		attr(res,"midpoint") <- midpoint
	}
	class(res) <- "dendrogram"
	return(res)
}

#' getMemberAttr
#'
#' fetch a given attribute from all leafs/members of a dendrogram
#'
#' @param tree		a dendrogram
#' @param attrName	the name of the attribute to be fetched
#' @return vector of attribute with one element for each leaf/member
#'
#' @author Fabian Mueller
#' @noRd
#' @examples
#' \donttest{
#' dhc <- as.dendrogram(hc <- hclust(dist(USArrests), "ave"))
#' dend <- dhc[[2]][[1]]
#' getMemberAttr(dend, "label")
#' }
getMemberAttr <- function(tree, attrName){
	Nmem <- attr(tree,"members")
	res <- dendrapply(tree, function(x){
		if (is.leaf(x)){
			return(attr(x, attrName))
		} else {
			return(NULL)
		}
	})
	return(unlist(res))
}


#' setMemberAttr
#'
#' set a given attribute to all leafs/members of a dendrogram
#'
#' @param tree		a dendrogram
#' @param attrName	the name of the attribute to be set
#' @param attrVals	a vector of values to be set for the leafs
#' @return a dendrogram in which the attributes for the leafs/members have been set to correspong values
#'
#' @author Fabian Mueller
#' @noRd
#' @examples
#' \donttest{
#' dhc <- as.dendrogram(hc <- hclust(dist(USArrests), "ave"))
#' dend <- dhc[[2]][[1]]
#' #assign some value to leafs
#' N <- attr(dend,"members")
#' res <- setMemberAttribute(dend, "someNum", sample(1:10,N,replace=TRUE))
#' str(unclass(res))
#' }
setMemberAttr <- function(tree, attrName, attrVals){
	Nmem <- attr(tree,"members")
	if (length(attrVals)==1){
		attrVals <- rep(attrVals,Nmem)
	}
	if (length(attrVals)!=Nmem){
		stop("Number of attribute values must be 1 or the number of tree members")
	}
	i <- 0
	res <- dendrapply(tree, function(x){
		if (is.leaf(x)){
			i <<- i + 1
			attr(x, attrName) <- attrVals[i]
		}
		return(x)
	})
	return(res)
}


#' dendrogramAttributeCombineRecursive
#'
#' recursively combine attributes using a defined combination function
#'
#' @param tree		a dendrogram
#' @param attrName	the name of the attribute to be set
#' @param combineFun a function definition that is used to combine the attributes for each level
#' @return a dendrogram in which the attributes for the leafs/members have been set to correspong values
#'
#' @details
#' the attribute values of all branches of a tree are recursively combined using the supplied function. The function itself takes as arguments the
#' attribute values from the daughter trees (not a vector of them)
#'
#' @author Fabian Mueller
#' @noRd
#' @examples
#' \donttest{
#' dhc <- as.dendrogram(hc <- hclust(dist(USArrests), "ave"))
#' dend <- dhc[[2]][[1]]
#' #assign some value to leafs
#' N <- attr(dend,"members")
#' res <- setMemberAttribute(dend, "someNum", sample(1:10,N,replace=TRUE))
#' rr <- dendrogramAttributeCombineRecursive(dend, "someNum")
#' str(unclass(rr))
#' rr <- dendrogramAttributeCombineRecursive(dend, "someNum", combineFun=sum)
#' str(unclass(rr))
#' }
dendrogramAttributeCombineRecursive <- function(tree, attrName, combineFun=c){
	if (is.leaf(tree)){
		return(tree)
	} else {
		subTreeList <- lapply(tree, FUN=function(x){
			dendrogramAttributeCombineRecursive(x, attrName, combineFun)
		})
		remergedTree <- tree
		#handle the case where there is just 1 subbranch of the tree
		if (length(subTreeList) > 1){
			# mp <- attr(tree,"midpoint") #somehow during merge, the midpoint gets lost. Save it and reset it later
			# remergedTree <- do.call("merge", c(subTreeList, list(height=attr(tree,"height"))))
			# attr(remergedTree,"midpoint") <- mp
			remergedTree <- dendrogramMergeSimple(subTreeList, height=attr(tree,"height"), midpoint=attr(tree,"midpoint"))
		}
		attrList <- lapply(remergedTree, FUN=function(x){attr(x, attrName)})
		attr(remergedTree, attrName) <- do.call(combineFun, attrList)
		return(remergedTree)
	}
}

