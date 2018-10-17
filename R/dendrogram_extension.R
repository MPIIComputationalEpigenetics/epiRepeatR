#' isLeaf
#'
#' helper function to determine if a dendrogram is a leaf, handling \code{NULL}
#' (in which case \code{FALSE} is returned)
#'
#' @param dend		a dendrogram
#' @return logical specifying whether the dendrogram is a leaf
#'
#' @author Fabian Mueller
#' @noRd
isLeaf <- function(dend){
	if (!is.null(attr(dend,"leaf"))) return(attr(dend,"leaf"))
	return(FALSE)
}

#' dendrogramMergeSimple
#'
#'  naively merge a list of dendrograms into a new dendrogram with each branch belonging to an element in the list by just setting the class and a couple of attributes
#'
#' @param dendroList	list of dendrograms
#' @param height		height (will be set as attribute) on which to merge the dendrograms
#' @param midpoint		midpoint attribute that will be set in the resulting dendrogram
#' @param label			label attribute that will be assigned to the resulting dendrogram
#' @return a dendrogram object with the dendrograms from the input list as branches
#'
#' @author Fabian Mueller
#' @noRd
dendrogramMergeSimple <- function(dendroList, height=NA, midpoint=NULL, label=NULL){
	res <- dendroList
	subTreeMems <- sapply(res,FUN=function(x){attr(x,"members")})
	attr(res,"members") <- sum(subTreeMems)
	attr(res,"height") <- height
	if (is.null(midpoint)){
		attr(res,"midpoint") <- (sum(subTreeMems)-1)/2
	} else {
		attr(res,"midpoint") <- midpoint
	}
	if (!is.null(label)) attr(res,"label") <- label
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

# dd <- list2dend(treeList)
# plot(dd)

#' createStump
#' 
#' create a stump (tree consisting of one leaf node)
#'
#' @param x		some label that gets turned into a stump
#'
#' @return a dendrogram for the stump
#'
#' @author Fabian Mueller
#' @noRd
createStump <- function(x){
	res <- x
	attr(res,"members") <- 1L
	attr(res,"label") <- x
	attr(res,"height") <- 0L
	attr(res,"leaf") <- TRUE
	class(res) <- "dendrogram"
	return(res)
}
#' createLinearTree
#' 
#' given a dendrogram, append a father node
#'
#' @param child      the dendrogram that becomes the child of the new father node
#' @param lab        label for the father node
#' @param height     height attribute for the father node
#'
#' @return a dendrogram with a new father node and the old dendrogram as its
#'         only child
#'
#' @author Fabian Mueller
#' @noRd
createLinearTree <- function(child, lab, h=(attr(child, "height")+1L)){
	res <- list(child)
	attr(res, "members")   <- attr(child, "members")
	attr(res, "height")    <- h
	attr(res, "label")     <- lab
	attr(res, "midpoint")  <- 0
	class(res) <- "dendrogram"
	return(res)
}
#' list2dend
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
		return(createStump(lab))
	} else if (length(x)==1){
		res <- createLinearTree(list2dend(x[[1]],names(x)[1]), lab)
	} else {
		subTrees <- lapply(names(x), FUN=function(ll){list2dend(x[[ll]], ll)})
		# stMembers <- vapply(subTrees, FUN=function(y){attr(y,"members")}, integer(1))
		stHeights <- vapply(subTrees, FUN=function(y){attr(y,"height")}, integer(1))

		res <- do.call("merge", unname(subTrees)) #careful: merge aparently modifies the midpoint attributes of subtrees again

		# attr(res,"members")  <- sum(stMembers)
		attr(res,"label")    <- lab
		attr(res,"height")   <- max(stHeights) + 1L
		attr(res,"midpoint") <- (attr(res,"members")-1)/2
	}
	return(res)
}
#' adjustAttr.midpoint
#'
#' recursively adjust the "midpoint" attribute
#' 
#' @param tree  the dendrogram to which the adjustment should be applied
#'
#' @return the modified dendrogram
#'
#' @author Fabian Mueller
#' @noRd
adjustAttr.midpoint <- function(tree){
	res <- dendrapply(tree, function(x){
		nMems <- attr(x,"members")
		attr(x, "midpoint") <- ifelse(nMems > 1, (nMems-1)/2, 0)
		return(x)
	})
	return(res)
}

#' makeBinary
#'
#' convert a dendrogram to a binary tree, i.e. each node is either a leaf or has exactly two children
#' (recursive function)
#' 
#' @param tree  the dendrogram to which the adjustment should be applied
#' @param epsilon a small number that will be used for pseudosplitting nodes with more than two children at different heights
#'
#' @return the modified dendrogram
#'
#' @author Fabian Mueller
#' @noRd
makeBinary <- function(tree, eps=1e-6){
	if (is.leaf(tree)) return(tree)
	if (length(tree)==2){
		return(dendrogramMergeSimple(list(makeBinary(tree[[1]], eps=eps), makeBinary(tree[[2]], eps=eps)), height=attr(tree,"height"), midpoint=attr(tree,"midpoint"), label=attr(tree,"label")))
	}
	if (length(tree)==1){
		#only 1 member: traverse down until the tree is either a leaf or not linear any more
		subTree <- tree[[1]] #dendrogramMergeSimple(tree[1],  height=attr(tree,"height"), midpoint=attr(tree,"midpoint"), label=NULL)
		attr(subTree, "height") <- attr(tree,"height") # Re-adjust the height
		attr(subTree, "midpoint") <- attr(tree,"midpoint") # Re-adjust the midpoint
		res <- makeBinary(subTree, eps=eps)
		return(res)
	}
	# more than 2 members
	tree.left <- makeBinary(tree[[1]])
	nMems.left <- attr(tree.left,"members")
	attr(tree.left, "midpoint") <- ifelse(nMems.left > 1, (nMems.left-1)/2, 0)

	tree.right <- makeBinary(dendrogramMergeSimple(tree[2:length(tree)], height=attr(tree,"height")-eps, label=attr(tree,"label")), eps=eps)
	nMems.right <- attr(tree.right,"members")
	attr(tree.right, "midpoint") <- ifelse(nMems.right > 1, (nMems.right-1)/2, 0)

	res <- dendrogramMergeSimple(list(tree.left, tree.right), height=attr(tree,"height"), midpoint=attr(tree,"midpoint"), label=attr(tree,"label"))
	return(res)
}
