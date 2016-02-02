library(devtools)
install_github("jtlovell/physGenomicsPVFinal")
library(physGenomicsPVFinal)
rm(list=ls())
data(potsGreenhouse)

info<-infoPots
counts<-countsPots
id.recovery<-info$id[info$Treatment=="recovery"]
sample.recovery<-info$Sample[info$id %in% id.recovery]
info.recovery<-info[info$Sample %in% sample.recovery, ]
counts.recovery<-counts[,sample.recovery]

wd<-which(info$Treatment!="recovery")
info.wetdry<-info[wd,]
counts.wetdry<-counts[,wd]
with(info.wetdry, table(Treatment, Time, Day))

#################################
# Part 4.1: Full model
#################################
stats<-pipeLIMMA(counts=counts.wetdry, info=info.wetdry, block=info$id, formula="~ Treatment * Time + Day")
voom.trtByMonth<-stats$voom[["E"]]
stats.fullmodel<-stats$simpleStats
stats.allests<-stats$stats
