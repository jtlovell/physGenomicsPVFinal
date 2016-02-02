rm(list=ls())
counts<-fread("/Users/John/Desktop/dropbox/Switchgrass_PlantPhys/R_package/switchgrassPhysGenomics/data/Eli/colorspace_GSM_counts.csv", stringsAsFactors=F)
# read in all the data for the shelter experiment
info<-read.csv("/Users/John/Desktop/dropbox/Switchgrass_PlantPhys/R_package/switchgrassPhysGenomics/data/Eli/colorspace_GSM_counts_X.csv", stringsAsFactors=F)
info$id<-sapply(info$Sample_title, function(x) strsplit(x, "_")[[1]][4])

info$id[info$id %in% c("8a","8b")]<-"8"
info$Day<-gsub(" ","", info$Day)
counts<-data.frame(counts)
rownames(counts)<-counts$gene
counts<-counts[,-1]
libs<-colSums(counts)

counts<-counts[,libs>1005000]
counts<-counts[(rowSums(counts)/ncol(counts))>=5,]

infoPots<-info[info$Sample %in% colnames(counts),]

countsPots<-counts[grepl("Pav",rownames(counts)),]

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
with(info.recovery, table(Treatment, Time, Day))
with(info.recovery, table(Treatment, Time, Day))

info.wetdry$id<-factor(info.wetdry$id, levels=unique(info.wetdry$id)[order(as.numeric(unique(info.wetdry$id)))])
info.recovery$id<-factor(info.recovery$id, levels=unique(info.recovery$id)[order(as.numeric(unique(info.recovery$id)))])

info.wetdry$Day<-factor(info.wetdry$Day,levels=c("13d","14d"))
info.recovery$Day<-factor(info.recovery$Day,levels=c("13d","14d"))

info.wetdry$Treatment<-factor(info.wetdry$Treatment,levels=c("drought","control"))
info.recovery$Treatment<-factor(info.recovery$Treatment,levels=c("drought","recovery"))

info.wetdry$Time<-factor(info.wetdry$Time,levels=c("5AM", "10AM", "12PM",  "2PM"))
info.recovery$Time<-factor(info.recovery$Time,levels=c("5AM", "10AM", "12PM",  "2PM"))

save(info.wetdry, info.recovery, counts.wetdry, counts.recovery, file="/Users/John/Desktop/dropbox/Switchgrass_PlantPhys/physGenomicsPVFinal/data/potsGreenhouse.RData")
