## Class for benhur object

setClass("BenHur", representation(jaccards = "list",
                                  size = "vector",
                                  iterations  = "vector",
                                  freq = "vector"),
         prototype = list(jaccards = list(),
           size = vector(), iterations = vector(),
           freq = vector()))

## Define a class for clusterComp

setClass("ClusterComp", representation(clusters = "vector",
                                       percent = "vector",
                                       freq = "vector",
                                       clusternum = "vector",
                                       iterations = "vector",
                                       method = "vector"),
         prototype = list(clusters = vector(), percent = vector(),
           freq = vector(), clusternum = vector(), iterations = vector(),
           method = vector()))

