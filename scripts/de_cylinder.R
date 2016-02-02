library(physGenomicsPVFinal)
rm(list=ls())

# this is the 2011 Pickle Cylinder data
data(cylinder)

# rename for convenience
counts<-countscyl
info<-infocyl

#################################
# Part 1: Full model
#################################
stats<-pipeLIMMA(counts=counts, info=info, block=info$pot, formula="~ trt * time")
voom.trtByMonth<-stats$voom[["E"]]
stats.fullmodel<-stats$simpleStats
stats.allests<-stats$stats

#################################
# Part 2: All pairwise contrasts
#################################
## Generate a design matrix that specifies all relavant contrasts
f<-info$time_trt
design <- model.matrix(~0+f)
contrast.matrix<-makeContrasts(fdawn_Wet-fdawn_Dry, fmidday_Wet-fmidday_Dry, fdusk_Wet-fdusk_Dry,
                               levels=design)
lim.contrasts<-anovaLIMMA(counts=counts, design=model.matrix(~0+f), block=info$pot, contrast.matrix=contrast.matrix)
#################################
# Part 3: Run PCA
#################################
pca<-voom2PCA(v=voom.trtByMonth, info=info, ids=info$id)
ggplot(pca, aes(x=PC1, y=PC2, col=time, shape=trt))+
  theme_bw()+geom_point(size=4)+scale_shape_manual(values=c(2,19))

#################################
# Part 4: Run model with MDWP as the predictor
#################################
stats<-pipeLIMMA(counts=counts, info=info, block=info$pot, formula="~ psi + time")
stats.fullmodel.mdwp<-stats$simpleStats
stats.allests.mdwp<-stats$stats

#################################
# Part 5: Run model with MDWP as the predictor, without controlling for location or year
#################################
stats<-pipeLIMMA(counts=counts, info=info, block=info$pot, formula="~ psi")
stats.fullmodel.mdwponly<-stats$simpleStats
stats.allests.mdwponly<-stats$stats

save(stats.fullmodel, stats.allests, lim.contrasts, pca, voom.trtByMonth,
     stats.fullmodel.mdwp, stats.allests.mdwp, stats.fullmodel.mdwponly, stats.fullmodel.mdwponly,
     file="/Users/John/Desktop/dropbox/Switchgrass_PlantPhys/stats_output/cylinder_allstats.RData")

