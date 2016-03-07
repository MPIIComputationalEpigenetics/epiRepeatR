#' treefork
#'
#' connect a root point to target points via a tree-like structure. Extension to the \code{diagram} package functionality
#'
#' @param from	the root point. A vector (or a one-row matrix) of one x and one y coordinate.
#' @param to	a matrix containing in its first column the x and in its second column the y coordinates for the target points
#' @param forkAt either "root" to create U-shaped strutures or "mid" to create forking (T-shaped) structures
#' @param path	split horizontally ("H") or vertically ("V")
#' @param lwd	line width
#' @param lty	line type
#' @param lcol  line color
#' @param ...	arguments passed down to dependent functions
#' @return a dendrogram in which the attributes for the leafs/members have been set to correspong values
#'
#' @details
#' modification of the slighly buggy treearrow function from package diagram
#'
#' @author Fabian Mueller
#' @noRd
#' @examples
#' \donttest{
#' openplotmat()
#' treefork(from=c(0,0.75), to=rbind(c(1,0.3), c(0.5,0.5), c(1,1)), forkAt="root", path="H")
#' openplotmat()
#' treefork(from=c(0.5,0.75), to=rbind(c(0,0.3), c(0.3,0.1), c(1,0)), forkAt="mid", path="V")
#' }
treefork <- function(from, to, forkAt="root", path="H", lwd=2, lty=1, lcol="black", ...){
	from <- matrix(from, ncol=2)
	to <- matrix(to, ncol=2)
	fork.x <- from[1,1]
	fork.hi.x <- from[1,1]
	fork.lo.x <- from[1,1]
	fork.y <- from[1,2]
	fork.hi.y <- from[1,2]
	fork.lo.y <- from[1,2]
	from.fork.xs <- to[,1]
	from.fork.ys <- to[,2]
	if (path=="H"){
		if (forkAt == "root"){
			fork.x <- from[1,1]
		} else if (forkAt == "mid"){
			fork.x <- mean(c(from[1,1], min(to[,1])))
		} else {
			stop("invalid forkAt argument")
		}
		
		fork.hi.x <- fork.x
		fork.lo.x <- fork.x
		fork.hi.y <- max(c(from[1,2], to[,2]))
		fork.lo.y <- min(c(from[1,2], to[,2]))
		from.fork.xs <- rep(fork.x, nrow(to))
		from.fork.ys <- to[,2]

	} else if (path=="V"){
		if (forkAt == "root"){
			fork.y <- from[1,2]
		} else if (forkAt == "mid"){
			fork.y <- mean(c(from[1,2], max(to[,2])))
		} else {
			stop("invalid forkAt argument")
		}
		fork.hi.y <- fork.y
		fork.lo.y <- fork.y
		fork.hi.x <- max(c(from[1,1], to[,1]))
		fork.lo.x <- min(c(from[1,1], to[,1]))
		from.fork.xs <- to[,1]
		from.fork.ys <- rep(fork.y, nrow(to))
	} else {
		stop("invalid path argument")
	}
	# straightarrow(from[1,], cbind(fork.x, fork.y), arr.side=0, ...)
	# straightarrow(cbind(fork.hi.x, fork.hi.y), cbind(fork.lo.x, fork.lo.y), arr.side=0, ...)
	# for (i in 1:length(from.fork.xs)){
	# 	straightarrow(cbind(from.fork.xs[i], from.fork.ys[i]), to[i,], arr.side=0, ...)
	# }
	x0 <- c(fork.hi.x, from.fork.xs)
	y0 <- c(fork.hi.y, from.fork.ys)
	x1 <- c(fork.lo.x, to[,1])
	y1 <- c(fork.lo.y, to[,2])
	if (forkAt != "root"){
		x0 <- c(from[1,1], x0)
		y0 <- c(from[1,2], y0)
		x1 <- c(fork.x, x1)
		y1 <- c(fork.y, y1)
	}
	segments(
		x0, y0, x1, y1,
		lwd=lwd, lty=lty, col=lcol
	)
}

