library(devtools)
install_github("jtlovell/physGenomicsPVFinal")
library(physGenomicsPVFinal)
data(temple2012_6treatments)

### for conviencence, rename info and counts
info<-info12
counts<-counts12

#################################
# Part 1: Full model with Treatment
#################################
## Run LIMMA, with subplot treated as repeated measures of the same plant.
### The only response is treatment. Sub_Block is the sub-block in the split plot design, containing two AP13
###    plants. The duplicate correlation and blocking in LIMMA should account for correlations between these plants
### We do not estimate intercept, so that the model only looks at the effects of Treatment, and the F-test reflects this
## Run the modeling pipeline
stats<-pipeLIMMA(counts=counts, info=info, block=info$Sub_Block, formula="~ Treatment")

### Extract voom-normalized counts
v<-stats$voom[["E"]]

### Extract statistics
### Here, the F statistics = "moderated F-statistics for testing all contrasts defined by the columns
###     of fit simultaneously equal to zero"
### See LIMMA::ebayes for descriptions of other statistics
stats.fullmodel<-stats$simpleStats

### Fstatistics for the whole model
stats.allests<-stats$stats

#################################
# Part 2: All pairwise contrasts
#################################
## Generate a design matrix that specifies all relavant contrasts
f<-info$Treatment
design <- model.matrix(~0+f)
contrast.matrix<-makeContrasts(f25th-flow, fmean-flow, fambient-flow, f75th-flow, fhigh-flow,
                               fmean-f25th, fambient-f25th, f75th-f25th, fhigh-f25th,
                               fambient-fmean, f75th-fmean, fhigh-fmean,
                               f75th-fambient, fhigh-fambient,
                               fhigh-f75th,
                               levels=design)
## run the LIMMA contrasts pipeline
## This returns a simple dataset with relevant statistics for each contrast, plus overall F-test.
lim.contrasts<-anovaLIMMA(counts=counts, design=design, block=info$Sub_Block, contrast.matrix=contrast.matrix)

#################################
# Part 3: Run PCA
#################################
pca<-voom2PCA(v=v, info=info, ids=info$ID)
ggplot(pca, aes(x=PC1, y=PC2, col=Treatment, alpha=Treatment))+
  theme_bw()+geom_point(size=4)

#################################
# Part 4: Subset to lines with physiology data
#################################
phys.lines<-info$ID[!is.na(info$order)]
counts<-counts[,phys.lines]
info<-info[info$ID %in% phys.lines,]
info$Treatment<-factor(info$Treatment, levels=c("low","mean","high"))

#################################
# Part 5: Run model with MDWP as the predictor
#################################
stats<-pipeLIMMA(counts=counts, info=info, block=info$Sub_Block, formula="~ MDWP")
stats.fullmodel.mdwp<-stats$simpleStats
stats.allests.mdwp<-stats$stats

#################################
# Part 6: Run model with Sampling Order as the predictor
#################################
stats<-pipeLIMMA(counts=counts, info=info, block=info$Sub_Block, formula="~ order")
stats.fullmodel.order<-stats$simpleStats
stats.allests.order<-stats$stats

#################################
# Part 7: Re-Run model for Treatment as the predictor
#################################
stats<-pipeLIMMA(counts=counts, info=info, block=info$Sub_Block, formula="~ Treatment")
stats.fullmodel.subtrt<-stats$simpleStats
stats.allests.subtrt<-stats$stats

save(stats.fullmodel, stats.allests, lim.contrasts, pca, v,
     stats.fullmodel.mdwp, stats.allests.mdwp, stats.fullmodel.order, stats.allests.order,stats.fullmodel.subtrt, stats.allests.subtrt,
     file="/Users/John/Desktop/dropbox/Switchgrass_PlantPhys/stats_output/tempe2012_allstats.RData")

