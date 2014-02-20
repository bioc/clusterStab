\name{ClusterComp-class}
\docType{class}
\alias{ClusterComp-class}
\alias{show,ClusterComp-method}

\title{Class "ClusterComp" a class for testing the stability of clusters
in microarray data}
\description{A specialized class representation used for testing the
  stability of clusters in microarray data.}
\section{Objects from the Class}{
Objects are usually created by a call to \code{clusterComp}, although
technically objects can be created by calls of the form
\code{new("ClusterComp", ...)}. However, the latter is probably not
worth doing.

}
\section{Slots}{
  \describe{
    \item{\code{clusters}:}{Object of class \code{"vector"} showing the
      cluster membership for each sample when using all the data. }
    \item{\code{percent}:}{Object of class \code{"vector"} containing
      the percentage of subsamples that resulted in the same class
      membership for all samples. }
    \item{\code{freq}:}{Object of class \code{"vector"} containing the
      subsampling percentage used. Defaults to 0.8. }
    \item{\code{clusternum}:}{Object of class \code{"vector"} containing
      the number of clusters tested.}
    \item{\code{iterations}:}{Object of class \code{"vector"} containing
    the number of iterations performed. Defaults to 100.}
    \item{\code{method}:}{Object of class \code{"vector"} containing the
      agglomerative method used. Options include "average", "centroid",
      "ward", "single", "mcquitty", or "median".}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "ClusterComp")}: Give a nice
      summary of results. }
  }
}
\references{  A. Ben-Hur, A. Elisseeff and I. Guyon. A stability based
  method for discovering structure in clustered data. Pacific Symposium
  on Biocomputing, 2002.
  Smolkin, M. and Ghosh, D. (2003).  Cluster stability scores for
  microarray data in cancer studies.
  BMC Bioinformatics 4, 36 - 42.
}
\author{James W. MacDonald <jmacdon@u.washington.edu> }
\keyword{classes}
