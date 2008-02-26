setMethod("clusterComp", "matrix", function(object, cl, seednum = NULL, B = 100,
                                            sub.frac = 0.8, method = "ave", adj.score = FALSE){

  do.clusterComp(object, cl, seednum = seednum, B = B, sub.frac = sub.frac, method = method,
                 adj.score = adj.score)
})

setMethod("clusterComp", "ExpressionSet", function(object, cl, seednum = NULL, B = 100,
                                             sub.frac = 0.8, method = "ave", adj.score = FALSE){
  
  object <- exprs(object)
  do.clusterComp(object, cl, seednum = seednum, B = B, sub.frac = sub.frac, method = method,
                 adj.score = adj.score)
})

## Show method

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
