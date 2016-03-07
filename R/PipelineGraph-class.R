setClass("PipelineGraph",
	slots = list(
		graph="igraph"
	),
	package = "epiRepeatR"
)
setMethod("initialize", "PipelineGraph",
	function(.Object){
		require(igraph)
		el <- rbind(
			c("seqReads", "bamExtract"),
			c("bamExtract", "repeatAlignment"),
			c("repeatAlignment", "methCalling"),
			c("repeatAlignment", "chipQuantification"),
			c("methCalling", "plotRepeatGroupTreesMeth"),
			c("methCalling", "plotRepeatMarkTree"),
			c("chipQuantification", "plotRepeatMarkTree")
		)
		.Object@graph <- graph_from_edgelist(el, directed = TRUE)
		.Object
	}
)
PipelineGraph <- function(){
	obj <- new("PipelineGraph")
	return(obj)
}

if (!isGeneric("getDependentSteps")) setGeneric("getDependentSteps", function(.Object, ...) standardGeneric("getDependentSteps"))
setMethod("getDependentSteps", signature(.Object="PipelineGraph"),
	function(.Object, step){
		if (!is.element(step, V(.Object@graph)$name)) stop("Invalid analysis step")
		res <- subcomponent(.Object@graph, step, mode="out")$name
		return(res)
	}
)
