######################
######################
# Part 1. Shelter data from June 2012, Temple, 6 treatments, AP13.
## Research Questions:  #1: What is the effect of time of sampling / MDWP?
##                      #2: What treatment contrasts give the strongest estimates?
##                      #3: What are the patterns of differential expression among treatments?
######################
######################

rm(list=ls())
setwd("/Users/JLovell/Dropbox/Switchgrass_PlantPhys/R_package/switchgrassPhysGenomics")
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

# 1.1: Read in the raw metadata, physiology data and counts.
# Metadata
layout<- read.csv("./data/shelter/TMPSH_Field Layout.csv", stringsAsFactors=F)
layout<-layout[layout$Genotype=="AP13",]
layout<-layout[,c(1,4,5,6)]; colnames(layout)[1]<-"PLOT"

basic<-read.delim("./data/shelter/TAGSeq_Shelter_AP12V2_all_counts_X.csv", stringsAsFactors=F)
basic<-basic[basic$GT=="AP13" ,c("ID","PLOT","Campaign","Treatment","Block")]
basic$Location<-ifelse(substr(basic$PLOT,1,1)==1, "WFC", "TMP")
basic$Year<-ifelse(grepl("2012",basic$Campaign),2012,ifelse(grepl("2013",basic$Campaign), 2013,2014))
basic$Month<-ifelse(grepl("June",basic$Campaign),"Jun","Aug")
table(basic$Location, basic$Year)
basic$Treatment<-tolower(basic$Treatment)

# Physiology
shelter.info<-read.csv("./data/shelter/shelter.info_sampleOrder.csv", stringsAsFactors=F)
shelter.info$Treatment<-tolower(shelter.info$Treatment)
shelter.info<-shelter.info[shelter.info$Genotype=="AP13",c("PLOT", "Year","Month","Location","Block","Treatment","Genotype","PDWP","MDWP","final.rna.order")]
shelter.info$Month<-ifelse(grepl("Jun",shelter.info$Month),"Jun","Aug")

# Counts
counts<-fread("./data/shelter/TAGSeq_Shelter_AP12V2_all_counts.csv")
setnames(counts, gsub("_merged","",names(counts)))
setnames(counts, gsub("_2$","",names(counts)))

# 1.2: Merge Phys and Metadata, cull to only 2012, temple, june, ap13
info<-merge(basic[,c("ID","PLOT", "Year","Month","Location","Treatment")],
            shelter.info[,c("PLOT", "Year", "Month", "PDWP", "MDWP","final.rna.order","Block")],
            by=c("PLOT", "Year","Month"), all.x=T)
samp.order.temple <- info[info$Location == "TMP"  & info$Year == 2012 & info$Treatment %in% c("high","low"),]
with(info, table(Treatment, Month, Year))
colnames(info)[which(colnames(info)=="final.rna.order")]<-"order"
info <- info[info$Year == 2012 & info$Month=="Jun",]
with(info, table(Treatment, Month, Year))

# 1.3: Cull the Counts data by IDs and Library Size
all.lines<-info$ID[info$ID %in% names(counts)]
counts<-counts[,c("gene",all.lines), with=F]
counts<-data.frame(counts)
norm.factors<-calcNormFactors(data.matrix(counts[,-1]))
lib.size<-colSums(counts[,-1])
# pdf("./stats_output/final_2012TMP/librarySize_plots.pdf")
plot(norm.factors, lib.size, main="comparison of library size and edgeR Normalization factors")
bad.lines<-names(counts)[-1][getOutliers(norm.factors)$iLeft]

counts<-counts[,-which(colnames(counts) %in% bad.lines)]
norm.factors<-calcNormFactors(data.matrix(counts[,-1]))
lib.size<-colSums(counts[,-1])
plot(norm.factors, lib.size, main="comparison of library size and edgeR Normalization factors",
     sub="following removal of libraries with very small normalization factors")
final.lines<-colnames(counts)[-1]

# 1.4: Cull genes by mean expression >5
means<-rowMeans(counts[,-1])
hist(log10(means+1), breaks=100, main="histogram of unfiltered gene mean expression \n (log10+1 scale)")
counts<-counts[means>5,]
means<-rowMeans(counts[,-1])
hist(log10(means+1), breaks=100, main="histogram of filter gene mean expression \n (log10+1 scale)")
rownames(counts)<-counts$gene
counts<-counts[,-1]
# dev.off()

# 1.5: Match information dataset with counts
info<-info[info$ID %in% final.lines,]
info<-info[match(final.lines,info$ID),]
cbind(info$ID, colnames(counts))

# write.csv(counts, "./stats_output/final_2012TMP/final.counts_temple2012.csv",row.names=F)
# write.csv(info, "./stats_output/final_2012TMP/final.info_temple2012.csv",row.names=F)

# 1.6: Re-set factors for expression analysis
info$Treatment<-factor(info$Treatment, levels=c("low","25th","mean", "ambient", "75th", "high"))
info$Treatment.num<-as.numeric(info$Treatment)
info<-merge(info[,-which(colnames(info)=="Block"),], layout, by="PLOT")
info12<-info
counts12<-counts
save(info12, counts12, file="/Users/JLovell/Dropbox/Switchgrass_PlantPhys/physGenomicsPVFinal/data/temple2012_6treatments.RData")
