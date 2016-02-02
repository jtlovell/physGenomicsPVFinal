library(devtools)
install_github("jtlovell/physGenomicsPVFinal")
library(physGenomicsPVFinal)
rm(list=ls())
data(potsGreenhouse)


#################################
# Part 4.1: Full model
#################################
stats<-pipeLIMMA(counts=counts.wetdry, info=info.wetdry, block=info$id, formula="~ Treatment * Time + Day")
voom.trtByMonth<-stats$voom[["E"]]
stats.fullmodel<-stats$simpleStats
stats.allests<-stats$stats
