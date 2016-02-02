library(physGenomicsPVFinal)
rm(list=ls())
data(tmpwfc20132014_2treatments)

### for conviencence, rename info and counts
info<-info1314
counts<-counts1314

#################################
# Part 2.1: Full model
#################################
stats<-pipeLIMMA(counts=counts, info=info, block=info$sb_unique, formula="~ Treatment * Location + Year + order")
v<-stats$voom[["E"]]
stats.fullmodel<-stats$simpleStats
stats.allests<-stats$stats

#################################
# Part 1.2: All pairwise contrasts
#################################
## Generate a design matrix that specifies all relavant contrasts

design <- model.matrix(~ 0+factor(Year_Location_Treatment.mat) + order , data=info)
colnames(design)[1:8]<-c("TMP_13_high","TMP_13_low","WFC_13_high","WFC_13_low",
                         "TMP_14_high","TMP_14_low","WFC_14_high","WFC_14_low")

contrast.matrix <- makeContrasts(TMP_13_high-TMP_13_low, WFC_13_high-WFC_13_low, TMP_14_high-TMP_14_low, WFC_14_high-WFC_14_low,
                                 levels=design)
lim.contrasts<-anovaLIMMA(counts=counts, design=design, block=info$Sub_Block, contrast.matrix=contrast.matrix)

#################################
# Part 3: Run PCA
#################################
pca<-voom2PCA(v=v, info=info, ids=info$ID)
ggplot(pca, aes(x=PC1, y=PC2, col=Location, shape=Treatment))+geom_point(size=4)+ facet_wrap(~Year, nrow=2)+
  theme_bw()+scale_shape_manual(values=c(2,19))

#################################
# Part 4: Run model with MDWP as the predictor
#################################
stats<-pipeLIMMA(counts=counts, info=info, block=info$sb_unique, formula="~ MDWP * Location + Year + order")
stats.fullmodel.mdwp<-stats$simpleStats
stats.allests.mdwp<-stats$stats

#################################
# Part 5: Run model with MDWP as the predictor, without controlling for location or year
#################################
stats<-pipeLIMMA(counts=counts, info=info, block=info$sb_unique, formula="~ MDWP + order")
stats.fullmodel.mdwponly<-stats$simpleStats
stats.allests.mdwponly<-stats$stats

save(stats.fullmodel, stats.allests, lim.contrasts, pca, v,
     stats.fullmodel.mdwp, stats.allests.mdwp, stats.fullmodel.mdwponly, stats.fullmodel.mdwponly,
     file="/Users/John/Desktop/dropbox/Switchgrass_PlantPhys/stats_output/shelter201314_allstats.RData")

