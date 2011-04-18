setMethod("benhur", "matrix", function(object, freq, upper, seednum = NULL, linkmeth = "average",
                                       distmeth = "euclidean", iterations = 100){

  do.benhur(object, freq, upper, seednum = seednum, linkmeth = linkmeth,
            iterations = iterations)
})

setMethod("benhur", "ExpressionSet", function(object, freq, upper, seednum = NULL,
                                              linkmeth = "average", distmeth = "euclidean",
                                              iterations = 100){

  object <- exprs(object)
  do.benhur(object, freq, upper, seednum = seednum, linkmeth = linkmeth,
            distmeth = distmeth, iterations = iterations)
})


## Show method for benhur


setMethod("show","BenHur",
          function(object){
            cat("A BenHur object\n")
            cat("produced using", object@iterations, "iterations and\n")
            cat("a ", object@freq * 100, "% subsampling frequency,\n", sep="")
            cat("testing for clusters <=", length(object@jaccards) + 1,"\n")
          })

## Define a generic for plotting jaccards


setMethod("hist", "BenHur",
          function(x, ...){
            ## Set up graphics device
            oldpar <- par(no.readonly = TRUE)
            par(mfrow=c(x@size,4))

            dat <- x@jaccards
            xlim <- range(hist(dat[[length(dat)]], 25, plot = FALSE)$breaks)
            for(i in seq(along=dat)){
              hist(dat[[i]], 25, freq = TRUE, main=paste("k = ", i + 1), ylab = "",
                   xlab = "Frequency", xlim = xlim)
            }
            par(oldpar)
          })

## Define a generic for plotting jaccardmatrix


setMethod("ecdf", "BenHur",
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
              legend(0, 1, legend = 2:(length(dat) + 1), lty = 1, col = 1:length(dat),
                     title = "Clusters")
              par(new = FALSE)

            }else{
              stop("There are too few samples to produce a cdf plot.\n",
                   "Only the histograms can be plotted.\n")

            }
          })
