rm(list=ls())
setwd("/Users/John/Desktop/dropbox/Switchgrass_PlantPhys/R_package/switchgrassPhysGenomics")
#####################
libs<-c("plyr")
invisible(lapply(libs, function(x) {cat(x,"..."); library(x, character.only=T, verbose=F, warn.conflicts=F,quietly=T)} ))


layout<- read.csv("./data/shelter/TMPSH_Field Layout.csv", stringsAsFactors=F)
layout<-layout[layout$Genotype=="AP13",]
layout<-layout[,c(1,4,5,6)]; colnames(layout)[1]<-"PLOT"
ltmp<-layout

layout<- read.csv("./data/shelter/WFCSH_Field Layout.csv", stringsAsFactors=F)
layout<-layout[layout$Genotype=="AP13",]
layout<-layout[,c(1,4,5,6)]; colnames(layout)[1]<-"PLOT"
lwfc<-layout
ls<-rbind(ltmp,lwfc)

# 1.1: Read in the raw metadata, physiology data and counts.
# Metadata
basic<-read.delim("./data/shelter/TAGSeq_Shelter_AP12V2_all_counts_X.csv", stringsAsFactors=F)
basic<-basic[basic$GT=="AP13" ,c("ID","PLOT","Campaign","Treatment","Block","Sub_Block")]
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
info<-merge(basic[,c("ID","PLOT", "Year","Month","Location","Treatment","Block","Sub_Block")],
            shelter.info[,c("PLOT", "Year", "Month", "PDWP", "MDWP","final.rna.order")],
            by=c("PLOT", "Year","Month"), all.x=T)
with(info, table(Treatment, Month, Year))
colnames(info)[which(colnames(info)=="final.rna.order")]<-"order"
info <- info[info$Year %in% c(2013,2014),]
with(info, table(Treatment, Location, Year))
info$ID<-gsub("[.]1","",info$ID)
info$ID<-gsub("[.]2","",info$ID)

info$sb_unique<-1
for(i in unique(info$PLOT)){
  sp<-ls[ls$PLOT==i,]
  if(nrow(sp)!=0){
    info$sb_unique[info$PLOT==i]<-sp$Sub_Block
  }
}

# 1.3: Cull the Counts data by IDs
all.lines<-info$ID
counts<-counts[,all.lines, with=F]
counts<-data.frame(counts)

means<-rowMeans(counts)
hist(log10(means+1), breaks=100, main="histogram of unfiltered gene mean expression \n (log10+1 scale)")
counts<-counts[means>5,]
means<-rowMeans(counts)
hist(log10(means+1), breaks=100, main="histogram of unfiltered gene mean expression \n (log10+1 scale)")

cbind(info$ID, colnames(counts))

# # 1.4: Cull genes by mean expression >5 and run PCA to cull outliers
# means<-rowMeans(counts[,-1])
# hist(log10(means+1), breaks=100, main="histogram of unfiltered gene mean expression \n (log10+1 scale)")
# counts<-counts[means>5,]
# means<-rowMeans(counts[,-1])
# hist(log10(means+1), breaks=100, main="histogram of unfiltered gene mean expression \n (log10+1 scale)")
#

# 1.6: Re-set factors for expression analysis
info$Treatment<-factor(info$Treatment, levels=c("high","mean","low"))
info$Year<-factor(info$Year, levels=c("2014","2013"))
info$Location<-factor(info$Location, levels=c("WFC","TMP"))
info$PLOT<-factor(info$PLOT)
info$LocTrtYear<-as.factor(with(info, paste(Location,Treatment,Year,sep="_")))
info$LocTrtYear.num<-as.numeric(info$LocTrtYear)
info<-ddply(info, .(Year, Month), mutate, order=rank(order))
info$sb_unique<-as.factor(info$Sub_Block)

save(info, counts, file="/Users/JLovell/Dropbox/Switchgrass_PlantPhys/physGenomicsPVFinal/data/tmpwfc20132014_3treatments.RData")
