\name{benhur}
\docType{methods}
\alias{benhur}
\alias{do.benhur}
\alias{benhur-methods}
\alias{benhur,matrix-method}
\alias{benhur,ExpressionSet-method}

\title{A Function to Estimate the Number of Clusters in Microarray Data }
\description{
 This function estimates the number of clusters in e.g., microarray data
 using an iterative process proposed by Asa Ben-Hur.
}
\usage{
\S4method{benhur}{ExpressionSet}(object, freq, upper, seednum = NULL,
linkmeth = "average", distmeth = "euclidean", iterations = 100)
\S4method{benhur}{matrix}(object, freq, upper, seednum = NULL, linkmeth
= "average", distmeth = "euclidean", iterations = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Either a matrix or \code{ExpressionSet} }
  \item{freq}{ The proportion of samples to use. This should be
    somewhere between 0.6 - 0.8 for best results.}
  \item{upper}{ The upper limit for number of clusters.}
  \item{seednum}{A value to pass to \code{set.seed}, which will allow
    for exact reproducibility at a later date.}
  \item{linkmeth}{Linkage method to pass to \code{hclust}. Valid values
    include "average", "centroid", "ward", "single", "mcquitty", or
    "median".}
  \item{distmeth}{The distance method to use. Valid values include
    "euclidean" and "pearson".}
  \item{iterations}{The number of iterations to use. The default of 100 is a
    reasonable number. }
}
\details{
  This function may be used to estimate the number of true clusters that
  exist in a set of microarray data. This estimate can be used to as
  input for \code{clusterComp} to estimate the stability of the clusters.

  The primary output from this function is a set of histograms that show
  for each cluster size how often similar clusters are formed from
  subsets of the data. As the number of clusters increases, the pairwise
  similarity of cluster membership will decrease. The basic idea is to
  choose the histogram corresponding to the largest number of clusters
  in which the majority of the data in the histogram is concentrated at
  or near 1.

  If overlay is set to \code{TRUE}, an additional CDF plot will be
  produced. This can be used in conjunction with the histograms to
  determine at which cluster number the data are no longer concentrated
  at or near 1.
}
\value{
  
  The output from this function is an object of class \code{benhur}. See
  the \code{benhur-class} man page for more information.
  
}
\references{ A. Ben-Hur, A. Elisseeff and I. Guyon. A stability based
  method for discovering structure in clustered data. Pacific Symposium
  on Biocomputing, 2002.
  Smolkin, M. and Ghosh, D. (2003).  Cluster stability scores for
  microarray data in cancer studies . BMC Bioinformatics 4, 36 - 42.
}
\author{Originally written by Mark Smolkin <marksmolkin@hotmail.com>
  further modifications by James W. MacDonald <jmacdon@med.umich.edu>}

\examples{
data(sample.ExpressionSet)
tmp <- benhur(sample.ExpressionSet, 0.7, 5)
hist(tmp)
ecdf(tmp)
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{ cluster }% __ONLY ONE__ keyword per line
