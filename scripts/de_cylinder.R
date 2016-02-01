library(physGenomicsPVFinal)
rm(list=ls())
data(cylinder)
# toss the column
counts<-counts[,-1]
counts<-counts[(rowSums(counts)/ncol(counts))>=5,]

#################################
# Part 2.1: Full model
#################################
stats<-pipeLIMMA(counts=counts, info=info, block=info$pot, formula="~ trt * time")
voom.trtByMonth<-stats$voom[["E"]]
stats.fullmodel<-stats$simpleStats
stats.allests<-stats$stats

#################################
# Part 1.2: All pairwise contrasts
#################################
## Generate a design matrix that specifies all relavant contrasts
f<-info$time_trt
design <- model.matrix(~0+f)
contrast.matrix<-makeContrasts(fdawn_Wet-fdawn_Dry, fmidday_Wet-fmidday_Dry, fdusk_Wet-fdusk_Dry,
                               levels=design)
lim.contrasts<-anovaLIMMA(counts=counts, design=model.matrix(~0+f), block=info$pot, contrast.matrix=contrast.matrix)
#################################
# Part 1.3: Run PCA
#################################
pca<-DESeq2PCA(counts=counts, info=info, formula="~ trt * time",
               factors2Plot=c("trt", "time"),
               factors2Extract=c("trt", "time","psi",  "RWC", "LDMC", "spad", "psi.delta", "Photo","Cond", "Ci","id","pot"))

#################################
# Part 1.4: Run model with MDWP as the predictor
#################################
stats<-pipeLIMMA(counts=counts, info=info, block=info$pot, formula="~ psi + time")
stats.fullmodel.mdwp<-stats$simpleStats
stats.allests.mdwp<-stats$stats

#################################
# Part 1.5: Run model with MDWP as the predictor, without controlling for location or year
#################################
stats<-pipeLIMMA(counts=counts, info=info, block=info$pot, formula="~ psi")
stats.fullmodel.mdwponly<-stats$simpleStats
stats.allests.mdwponly<-stats$stats

save(stats.fullmodel, stats.allests, lim.contrasts, pca, voom.trtByMonth,
     stats.fullmodel.mdwp, stats.allests.mdwp, stats.fullmodel.mdwponly, stats.fullmodel.mdwponly,
     file="/Users/John/Desktop/dropbox/Switchgrass_PlantPhys/stats_output/cylinder_allstats.RData")

