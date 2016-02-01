rm(list=ls())
setwd("/Users/John/Desktop/dropbox/Switchgrass_PlantPhys/R_package/switchgrassPhysGenomics")
opar<-par()
#####################
libs<-c("DESeq2", "limma", "edgeR","qvalue", # differential expression analyses
        "plyr", "tidyr", "reshape", "reshape2", "data.table", "qdap","gdata", # data manipulation
        "lme4", "MASS", "vegan","lmerTest", "car","extremevalues", # stats
        "ggplot2", "VennDiagram","venneuler","Heatplus","dendroextras")  # plotting
invisible(lapply(libs, function(x) {cat(x,"..."); library(x, character.only=T, verbose=F, warn.conflicts=F,quietly=T)} ))


source("./R/pipeLIMMA.R")
source("./R/pipeLIMMA.voomfit.R")
source("./R/DESeq2PCA.R")
source("./R/counts2Venn.R")
source("./R/getContrastStats.R")
source("./R/getEbayesStats.R")
source("./R/limmaStats.R")
source("./R/anovaLIMMA.R")
source("./R/voom2MeanHeatMaps.R")

#####################
#####################
# Part #1: Get information read in and screen for outliers

#####################
# 1.1: read in the counts and informations
# ****Note**** info MUST contain a column called "ID" that exactly matches the column names of counts,
# ****Note**** with the exception that the first column of counts MUST be called "gene", this is tossed in step 1.4
info<-read.csv("./data/all_final/shelter_tmp.wfc.2013.2014_info.csv")
counts<-read.csv("./data/all_final/shelter_tmp.wfc.2013.2014_counts.csv")

#####################
# 1.2: set the levels of all data to be analyzed (that are not continuous)
info$Year<-factor(info$Year, levels=c("2013","2014"))
info$Location<-factor(info$Location, levels=c("TMP","WFC"))
info$Treatment<-factor(info$Treatment, levels=c("high","low"))
info$Block<-as.factor(info$Block)

#####################
# 1.3: generate multi-level factors for use with contrasts
info$Year_Location_Treatment<-as.factor(with(info, paste(Year,Location,Treatment, sep="_")))
info$Year_Location_Treatment.mat<-as.numeric(info$Year_Location_Treatment)
#####################
# 1.4: cull counts so that they have exactly the same dimensions
row.names(counts)<-counts$gene
counts<-counts[,-1]

#####################
# 1.5 PCA to define outliers and toss if necessary
pca<-DESeq2PCA(counts=counts, info=info, formula="~ Treatment * Location * Year + order",
               factors2Plot=c("Treatment", "Location", "Year"),
               factors2Extract=c("Treatment", "Location", "Year","order","MDWP", "PDWP","ID"))
bad.lines<-as.character(pca$ID[pca$PC2 <= (-40)])

#####################
# 1.6: Cull counts and info so that they represent only the good lines
info<-info[!info$ID %in% bad.lines,]
counts<-counts[,-which(colnames(counts) %in% bad.lines)]

#####################
# 1.7: check that everything lines up
cbind(as.character(info$ID), colnames(counts))

layout<- read.csv("./data/shelter/TMPSH_Field Layout.csv", stringsAsFactors=F)
layout<-layout[layout$Genotype=="AP13",]
layout<-layout[,c(1,4,5,6)]; colnames(layout)[1]<-"PLOT"
ltmp<-layout

layout<- read.csv("./data/shelter/WFCSH_Field Layout.csv", stringsAsFactors=F)
layout<-layout[layout$Genotype=="AP13",]
layout<-layout[,c(1,4,5,6)]; colnames(layout)[1]<-"PLOT"
lwfc<-layout
ls<-rbind(ltmp,lwfc)

info$sb_unique<-1
for(i in unique(info$PLOT)){
  sp<-ls[ls$PLOT==i,]
  if(nrow(sp)!=0){
    info$sb_unique[info$PLOT==i]<-sp$Sub_Block
  }
}
info1314<-info
counts1314<-counts
save(info1314, counts1314, file="/Users/John/Desktop/dropbox/Switchgrass_PlantPhys/physGenomicsPVFinal/data/tmpwfc20132014_3treatments.RData")
