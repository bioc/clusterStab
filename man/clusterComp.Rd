\name{clusterComp}
\docType{methods}
\alias{clusterComp}
\alias{do.clusterComp}
\alias{clusterComp-methods}
\alias{clusterComp,matrix-method}
\alias{clusterComp,exprSet-method}


\title{ Estimate Microarray Cluster Stability}
\description{
 This function estimates the stability of clustering solutions using
 microarray data. Currently only agglomerative hierarchical clustering
 is supported.
}
\usage{
\S4method{clusterComp}{exprSet}(object, cl, seednum = NULL, B = 100,
sub.frac = 0.8, method = "ave", adj.score = FALSE)
\S4method{clusterComp}{matrix}(object, cl, seednum = NULL, B = 100,
sub.frac = 0.8, method = "ave", adj.score = FALSE)
}
\arguments{
  \item{object}{Either a matrix or \code{exprSet} }
  \item{cl}{ The number of clusters. This may be estimated using \code{benhur}}
  \item{seednum}{A value to pass to \code{set.seed}, which will allow
    for exact reproducibility at a later date.}
  \item{B}{ The number of permutations.}
  \item{sub.frac}{The proportion of genes to use in each
    subsample. This value should be in the range of 0.75 - 0.85 for best
    results}
  \item{method}{ The linkage method to pass to \code{hclust}.  Valid values
    include "average", "centroid", "ward", "single", "mcquitty", or
    "median".}
  \item{adj.score}{Boolean. Should the stability scores be adjusted for
    cluster size? Defaults to \code{FALSE}.}
}
\details{
  This function estimates the stability of a clustering solution by
  repeatedly subsampling the data and comparing the cluster membership
  of the subsamples to the original clusters.
}
\value{
  
  The output from this function is an object of class \code{clusterComp}. See
  the \code{clusterComp-class} man page for more information.
  
}
\references{ A. Ben-Hur, A. Elisseeff and I. Guyon. A stability based
  method for discovering structure in clustered data. Pacific Symposium
  on Biocomputing, 2002.
  Smolkin, M. and Ghosh, D. (2003).  Cluster stability scores
  for microarray data in cancer studies . BMC Bioinformatics 4, 36 - 42. }
\author{ James W. MacDonald <jmacdon@med.umich.edu>}
\examples{
data(sample.exprSet)
clusterComp(sample.exprSet, 3)
}
\keyword{cluster }% at least one, from doc/KEYWORDS
