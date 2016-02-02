library(physGenomicsPVFinal)
rm(list=ls())
data(potsGreenhouse)

#################################
# Part 4.1: Full model
#################################
stats<-pipeLIMMA(counts=countsPots, info=info, block=info$pot, formula="~ trt * time")
voom.trtByMonth<-stats$voom[["E"]]
stats.fullmodel<-stats$simpleStats
stats.allests<-stats$stats
