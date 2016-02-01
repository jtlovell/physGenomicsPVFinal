library(physGenomicsPVFinal)
rm(list=ls())
data(tmpwfc20132014_3treatments)

### for conviencence, rename info and counts
info<-info1314
counts<-counts1314

#################################
# Part 2.1: Full model
#################################
stats<-pipeLIMMA(counts=counts, info=info, block=info$sb_unique, formula="~ Treatment * Location + Year + order")
voom.trtByMonth<-stats$voom[["E"]]
stats.fullmodel<-stats$simpleStats
stats.allests<-stats$stats

#################################
# Part 1.2: All pairwise contrasts
#################################
## Generate a design matrix that specifies all relavant contrasts
design <- model.matrix(~ 0+factor(Year_Location_Treatment.mat) + order , data=info)
colnames(design)[1:8]<-levels(info$Year_Location_Treatment)
colnames(design)[1:8]<-multigsub(c("2013","2014"),c("a","b"), colnames(design)[1:8])
contrast.matrix <- makeContrasts(a_TMP_high-a_TMP_low,  b_TMP_high-b_TMP_low,
                                 a_WFC_high-a_WFC_low, b_WFC_high-b_WFC_low,
                                 levels=design)
lim.contrasts<-anovaLIMMA(counts=counts, design=design, block=info$Sub_Block, contrast.matrix=contrast.matrix)

#################################
# Part 1.3: Run PCA
#################################
pca<-DESeq2PCA(counts=counts, info=info, formula="~ Treatment * Location + Year + order",
               factors2Plot=c("Treatment", "Location", "Year"),
               factors2Extract=c("Treatment", "Location", "Year","order","MDWP", "PDWP","ID"))

#################################
# Part 1.4: Run model with MDWP as the predictor
#################################
stats<-pipeLIMMA(counts=counts, info=info, block=info$Sub_Block, formula="~ MDWP * Location + Year + order")
stats.fullmodel.mdwp<-stats$simpleStats
stats.allests.mdwp<-stats$stats

#################################
# Part 1.5: Run model with MDWP as the predictor, without controlling for location or year
#################################
stats<-pipeLIMMA(counts=counts, info=info, block=info$Sub_Block, formula="~ MDWP + order")
stats.fullmodel.mdwponly<-stats$simpleStats
stats.allests.mdwponly<-stats$stats


save(stats.fullmodel, stats.allests, lim.contrasts, pca, voom.trtByMonth,
     stats.fullmodel.mdwp, stats.allests.mdwp, stats.fullmodel.mdwponly, stats.fullmodel.mdwponly,
     file="/Users/John/Desktop/dropbox/Switchgrass_PlantPhys/stats_output/shelter201314_allstats.RData")

