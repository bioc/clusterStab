\name{BenHur-class}
\docType{class}
\alias{BenHur-class}
\alias{ecdf,BenHur-method}
\alias{hist,BenHur-method}
\alias{show,BenHur-method}

\title{Class "BenHur", a class for estimating clusters in microarray
  data, and methods for visualizing them.}
\description{A specialized class representation used for estimating
  clusters in microarray data.}
\section{Objects from the Class}{
Objects are usually created by a call to \code{benhur}, although
technically a new object can also be created by a call to
\code{new("BenHur",...)}. However, this second method is usually not
worth the work required.
}
\section{Slots}{
  \describe{
    \item{\code{jaccards}:}{Object of class \code{"list"}, containing the
    jaccard vectors; these indicate the proportion of pairwise
    similarity between clusters formed from subsets of the data.}
    \item{\code{size}:}{Object of class \code{"vector"}, only used for plotting.}
    \item{\code{iterations}:}{Object of class \code{"vector"},
      containing the number of iterations. Defaults to 100.}
    \item{\code{freq}:}{Object of class \code{"vector"}, containing the
      proportion of the data used for subsampling. }
  }
}
\section{Methods}{
  \describe{
    \item{ecdf}{\code{signature(x = "BenHur")}: Plot an empirical
      CDF. This can be used to help determine the number of clusters in
      the data. The most likely (e.g., most stable number) of clusters
      will have a CDF that is concentrated at or near one. See vignette
      for more information.}
    \item{hist}{\code{signature(x = "BenHur")}: Plot histograms for all
      clusters tested. The most likely (e.g., most stable number) of
      clusters will have a histogram in which the data are clustered at
      or near one. See vignette for more information.}
    \item{show}{\code{signature(object = "BenHur")}: Gives a nice
      summary. }
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
