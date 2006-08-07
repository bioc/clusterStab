############################################################3
##
## James W. MacDonald - 2005
##
##  clusterStab - Functions to estimate sample clusters and to
##                determine the stability of clusters in microarray data
##
## Modifications:
##
## 1-5-05 Cleaned up the original code and optimized for speed.
## 3-22-05 Added  S4 classes benhur and clusterComp, as well as methods for them.
## 2-28-05 Improved computation for Cone and Ctwo
## 6-13-05 Modified do.clusterComp to assess individual cluster stability
##         instead of overall cluster stability as before. In addition, added
##         the adjusted score calculation to account for cluster size.
##
####################################################################






###############################################################################################
##
## benhur - a function to estimate the number of sample clusters for microarray data, based on the methods in
## A. Ben-Hur, A. Elisseeff and I. Guyon. A stability based method for discovering structure in clustered data.
## Pacific Symposium on Biocomputing, 2002.
##
## Originally written by Mark Smolkin <marksmolkin@hotmail.com> further modifications by James W. MacDonald
## <jmacdon@med.umich.edu>
##
## Inputs:
##
## object - a matrix of data where columns are samples and rows are genes. These data should already be normalized
##
## freq - the proportion of samples to use. This value should exceed 0.5 (and of course be less than 1).
##
## upper - the upper limit for the number of clusters.
##
## seednum - an integer used to set the seed. This allows for an exact reproduction of the analysis at a later date.
##
## linkmeth - the linkage method to use for clustering. Valid values include "average", "centroid", "ward", "single",
##            "mcquitty", or "median".
##
## iterations - the number of subsamplings to perform.
##
##
## Output:
##
## An object of class benhur.
##
##
################################################################################################

do.benhur <- function(object, freq , upper, seednum = NULL, linkmeth = "average",  iterations = 100){
  
  if(!is.null(seednum))
    set.seed(seednum)
  require(stats, quietly = TRUE)
  
  ## Microarray data is usually samples in columns, genes in rows, so we must transpose
  dat <- t(object)
  size <- ceiling((upper-1)/4)
  samples <- dim(dat)[1]
  n <- dim(dat)[2]
  
  jaccardmatrix <- NULL
  
  
  ## Make a list to hold jaccard matrices
  
  jaccards <- vector("list", length = (upper - 1))
  
  
  for(m in 1:iterations){
    sampleone <- sample(1:samples, freq*samples, replace = FALSE)
    sampletwo <- sample(1:samples, freq*samples, replace = FALSE)
    
    subone <- dat[sampleone,]
    subtwo <- dat[sampletwo,]
    
    
    
    clustone <- hclust(dist(subone), method = linkmeth)
    clusttwo <- hclust(dist(subtwo), method = linkmeth)
    
    
    ##Checking for intersection
    
    inboth <- sampleone[sampleone %in% sampletwo]
    
    
    for(k in 2:upper){
      
      ## Cut the trees
      
      cutone <- cutree(clustone, k)
      cuttwo <- cutree(clusttwo, k)
      
      ## Cluster membership for intersecting samples
      
      ord <- 1:length(sampleone)
      
      placeone <- cutone[ord[match(inboth, sampleone)]]
      placetwo <- cuttwo[ord[match(inboth, sampletwo)]]
      
      
      ## This code forms the 'C' matrices s.t. Cij = 1 if xi and xj belong to the same cluster
      ## and i and j are not equal. Cij = 0 otherwise
      
      m.dim <- length(placeone)
      Cone <- matrix(placeone, nr = m.dim, nc = m.dim) == matrix(placeone, nr = m.dim, nc = m.dim, byrow = TRUE)
      Ctwo <- matrix(placetwo, nr = m.dim, nc = m.dim) == matrix(placetwo, nr = m.dim, nc = m.dim, byrow = TRUE)
      diag(Cone) <- diag(Ctwo) <- 0
      
      jaccard <- sum(Cone * Ctwo)/(sum(Cone) + sum(Ctwo) - sum(Cone * Ctwo))
      jaccards[[(k-1)]] <- c(jaccards[[(k-1)]], jaccard)
    }
  }
  
  ## There is a chance that some of the jaccard matrices consist of all ones, which will result
  ## in a histogram that is just a big box. To account for this possibility, we adjust the first value by 0.00001
  
  jac.test <- lapply(lapply(jaccards, function(x) table(x)), function(y) length(y)) == 1
  if(any(jac.test)){
    fix.jacs <- which(jac.test)
    for(i in fix.jacs)
      jaccards[[i]][[1]] <- jaccards[[i]][[1]] - 0.0001
  }

  ans <- new("benhur", jaccards = jaccards, size = size, iterations = iterations,
             freq = freq)
  ans
}
  
   
 

########################################################################################
##
## clusterComp - A function to compute stability of clusters based on subsampling as described
##  in  Smolkin, M. and Ghosh, D. (2003).  Cluster stability scores for microarray data in
##  cancer studies . BMC Bioinformatics 4, 36 - 42.
##
## Copyright James W. MacDonald, 2005, all rights reserved
## 
##
## Inputs:
##
## object - a matrix of expression data where columns are samples and rows are genes
##
## cl - the number of clusters
##
## seednum - an integer value used to set the random seed.http://svnbook.red-bean.com/svnbook-1.1/svn-book.html#svn-ch-9
##
## B - the number of iterations to perform
##
## sub.frac - the fraction of genes to use in each subsampling. Recommended values are between 0.75 and 0.85
##
## method - the agglomeration method to use. Choices include "ward", "single", "average", "complete", "mcquitty",
##          "median", or "centroid". See ?hclust for more information.
##
## Output:
##
## A list containing:
##
## Clusters - a vector identifying the cluster membership of each sample
##
## Percent - the percentage of subsampled clusters that were identical to the original cluster
##
## Clusternumber - the number of clusters
##
## Genesused - the fraction of genes used for each subsampling
##
## 
##
#############################################################################################

do.clusterComp <- function(object, cl, seednum = NULL, B = 100, sub.frac = 0.8, method = "ave",
                           adj.score = FALSE){
  
  if(!is.null(seednum))
    set.seed(seednum)

  methods <- c("ward", "single", "complete", "average", "mcquitty", 
               "median", "centroid")

  ## Make sure cl is an integer
  
  if(length(cl) != 1 || !is.numeric(cl))
    stop("The cl argument should only be a single number!")
  
  ## Get original cluster memberships
  orig.clust <- hclust(dist(t(object)), method = method)
  orig.memb <- cutree(orig.clust, cl)
  memb.list <- vector("list", length=cl)
  for(i in 1:cl)
    memb.list[[i]] <- which(orig.memb == i)
  
  ## Subsample and output cluster memberships in a matrix
  
  len <- dim(object)[1]
  count <- vector("list", cl)
  cnt <- 0
  while(cnt < B){
    index <- sample(1:len, sub.frac * len, replace = FALSE)
    tmp <- hclust(dist(t(object[index,])), method = method)
    memb <- cutree(tmp, cl)
    for(i in 1:cl){
      tst <- all(memb[memb.list[[i]]] == i)
      count[[i]] <- c(count[[i]], tst)
    }
    cnt <- cnt + 1
  }
  if(adj.score){
    sizes <- lapply(memb.list, function(x) 1/(log(length(x)) + 1))
    sc.adj <- function(x, y) ((sum(x)/B)^y)*100
    percent <- mapply(sc.adj, count, sizes)
    }else{
      percent <- sapply(count, function(x) sum(x)*100/B)
    }
  ans <- new("clusterComp", clusters = orig.memb, percent = round(percent,2),
             freq = sub.frac, clusternum = cl, iterations = B,
             method = match.arg(method, methods))
  ans
}

#####################################################################################
##
## S4 methods for the above functions
##
#####################################################################################

if(is.null(getGeneric("benhur")))
  setGeneric("benhur", function(object, ...) standardGeneric("benhur"))


setMethod("benhur", "matrix", function(object, freq, upper, seednum = NULL, linkmeth = "average",
                                       iterations = 100){
  
  do.benhur(object, freq, upper, seednum = seednum, linkmeth = linkmeth,
            iterations = iterations)
})

setMethod("benhur", "exprSet", function(object, freq, upper, seednum = NULL, linkmeth = "average",
                                        iterations = 100){
  
  object <- object@exprs
  do.benhur(object, freq, upper, seednum = seednum, linkmeth = linkmeth,
            iterations = iterations)
})

setMethod("benhur", "ExpressionSet", function(object, freq, upper, seednum = NULL, linkmeth = "average",
                                        iterations = 100){
  
  object <- exprs(object)
  do.benhur(object, freq, upper, seednum = seednum, linkmeth = linkmeth,
            iterations = iterations)
})


if(is.null(getGeneric("clusterComp")))
  setGeneric("clusterComp", function(object, ...) standardGeneric("clusterComp"))

setMethod("clusterComp", "matrix", function(object, cl, seednum = NULL, B = 100,
                                            sub.frac = 0.8, method = "ave", adj.score = FALSE){

  do.clusterComp(object, cl, seednum = seednum, B = B, sub.frac = sub.frac, method = method,
                 adj.score = adj.score)
})

setMethod("clusterComp", "exprSet", function(object, cl, seednum = NULL, B = 100,
                                             sub.frac = 0.8, method = "ave", adj.score = FALSE){
  
  object <- object@exprs
  do.clusterComp(object, cl, seednum = seednum, B = B, sub.frac = sub.frac, method = method,
                 adj.score = adj.score)
})

setMethod("clusterComp", "ExpressionSet", function(object, cl, seednum = NULL, B = 100,
                                             sub.frac = 0.8, method = "ave", adj.score = FALSE){
  
  object <- exprs(object)
  do.clusterComp(object, cl, seednum = seednum, B = B, sub.frac = sub.frac, method = method,
                 adj.score = adj.score)
})

## Class for benhur object

setClass("benhur", representation(jaccards = "list",
                                  size = "vector",
                                  iterations  = "vector",
                                  freq = "vector"),
         prototype = list(jaccards = list(),
           size = vector(), iterations = vector(),
           freq = vector()))

## Show method for benhur

if(is.null(getGeneric("show")))
  setGeneric("show")

setMethod("show","benhur",
          function(object){
            cat("A benhur object\n")
            cat("produced using", object@iterations, "iterations and\n")
            cat("a ", object@freq * 100, "% subsampling frequency,\n", sep="")
            cat("testing for clusters <=", length(object@jaccards) + 1,"\n")
          })

## Define a generic for plotting jaccards

if(is.null(getGeneric("hist")))
  setGeneric("hist", function(x, ...) standardGeneric("hist"))

setMethod("hist", "benhur",
          function(x, ...){
            ## Set up graphics device
            oldpar <- par(no.readonly = TRUE)
            par(mfrow=c(x@size,4))
            
            dat <- x@jaccards
            xlim <- range(hist(dat[[length(dat)]], 25, freq = TRUE, plot = FALSE)$breaks)
            for(i in seq(along=dat)){
              hist(dat[[i]], 25, freq = TRUE, main=paste("k = ", i + 1), ylab = "", xlab = "Frequency", xlim = xlim)
            }
            par(oldpar)
          })

## Define a generic for plotting jaccardmatrix

if(is.null(getGeneric("ecdf")))
  setGeneric("ecdf", function(x) standardGeneric("ecdf"))

setMethod("ecdf", "benhur",
          function(x){
            seqvec <- seq(1/x@iterations, 1, by = 1/x@iterations)
            dat <- x@jaccards
            ## Make sure all jaccards are of the same length
            if(all(x@iterations == sapply(dat, length))){
              for (i in 1:length(dat)){
                plot(sort(dat[[i]]), seqvec, type = "s", cex = 0.1, xlim = c(0, 1), ylim = c(0, 1),
                     ylab = "cumulative", xlab = "similarity", col = i)
                par(new = TRUE)
              }
              legend(0, 1, legend = 2:(length(dat) + 1), lty = 1, col = 1:length(dat), title = "Clusters")
              par(new = FALSE)
              
            }else{
              stop("There are too few samples to produce a cdf plot.\n Only the histograms can be plotted.\n")
              
            }
          })

## Define a class for clusterComp

setClass("clusterComp", representation(clusters = "vector",
                                       percent = "vector",
                                       freq = "vector",
                                       clusternum = "vector",
                                       iterations = "vector",
                                       method = "vector"),
         prototype = list(clusters = vector(), percent = vector(),
           freq = vector(), clusternum = vector(), iterations = vector(),
           method = vector()))

## Show method

if(is.null(getGeneric("show")))
  setGeneric("show")

setMethod("show","clusterComp",
          function(object){
            stab.out <- matrix(paste(round(object@percent, 0), "%", sep=""), nr = 1)
            dimnames(stab.out) <- list("Cluster stability: ", 1:length(object@percent))
            cat("Results from running clusterComp:\n")
            print(stab.out, quote = FALSE)
            cat("Iterations:", object@iterations, "\n")
            cat("Subsampling frequency: ", object@freq * 100, "%\n", sep="")
            cat("Agglomeration method:", object@method, "\n")
            cat("Original cluster membership:\n", object@clusters,"\n")
          })


    
