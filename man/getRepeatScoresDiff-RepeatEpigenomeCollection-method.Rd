% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RepeatEpigenomeCollection-class.R
\docType{methods}
\name{getRepeatScoresDiff,RepeatEpigenomeCollection-method}
\alias{getRepeatScoresDiff,RepeatEpigenomeCollection-method}
\alias{getRepeatScoresDiff}
\title{getRepeatScoresDiff-methods}
\usage{
\S4method{getRepeatScoresDiff}{RepeatEpigenomeCollection}(.Object, mark,
  compInfo)
}
\arguments{
\item{.Object}{\code{\linkS4class{RepeatEpigenomeCollection}} object}

\item{mark}{Epigenetic mark for which the score information is retrieved}

\item{compInfo}{\code{S3} object of class \code{comparisonInfo} as returned by \code{\link{getComparisonInfo,RepeatEpigenomeCollection-method}}
or \code{\link{getComparisons,RepeatEpigenomeCollection-method}}}
}
\value{
A matrix with differential statistics for each repeat subfamily.
}
\description{
Retrieve values of differential epigenetic quantifications for each RE in each sample
}
\author{
Fabian Mueller
}
