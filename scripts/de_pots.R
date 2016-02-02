library(devtools)
install_github("jtlovell/physGenomicsPVFinal")
library(physGenomicsPVFinal)
rm(list=ls())
data(potsGreenhouse)

info<-info.wetdry
counts<-counts.wetdry


library(qdap)
info$Treatment<-multigsub(c("control","drought","14d","13d","5AM","10AM","12PM","2PM"),
                          c("w","d","d1","d2","t0","t5","t7","t9"),
                          as.character(info$Treatment))
info$Time<-multigsub(c("control","drought","14d","13d","5AM","10AM","12PM","2PM"),
                          c("w","d","d1","d2","t0","t5","t7","t9"),
                          as.character(info$Time))
info$Day<-multigsub(c("control","drought","14d","13d","5AM","10AM","12PM","2PM"),
                          c("w","d","d1","d2","t0","t5","t7","t9"),
                          as.character(info$Day))
info$time_day_try<-as.factor(with(info, paste(Treatment,Day,Time, sep="_")))
info$Treatment<-factor(info$Treatment, levels=c("w","d"))
info$Time<-factor(info$Time, levels=c("t0", "t5", "t7", "t9"))
info$Day<-factor(info$Day, levels=c("d1","d2"))
with(info, table(Time, Day, Treatment))
#################################
# Part 1: Full model
#################################
stats<-pipeLIMMA(counts=counts, info=info, block=info$id, formula="~ Treatment * Time + Day")
v<-stats$voom[["E"]]
stats.fullmodel<-stats$simpleStats
stats.allests<-stats$stats

#################################
# Part 2: Run PCA
#################################
pca<-voom2PCA(v=v, info=info, ids=info$id)
library(ggplot2)
ggplot(pca, aes(x=PC1, y=PC2, col=Day, shape=Treatment))+theme_bw() +facet_wrap(~Time)+
  geom_point(size=4)+scale_shape_manual(values=c(2,19))

#################################
# Part 3: All pairwise contrasts
#################################
## Generate a design matrix that specifies all relavant contrasts
f<-info$time_day_try
design <- model.matrix(~0+f)
contrast.matrix<-makeContrasts(fw_d1_t0-fd_d1_t0, fw_d1_t5-fd_d1_t5, fw_d1_t7-fd_d1_t7,fw_d1_t9-fd_d1_t9,
                               fw_d2_t0- fd_d2_t0, fw_d2_t9- fd_d2_t9,
                               levels=design)
lim.contrasts<-anovaLIMMA(counts=counts, design=design, block=info$id, contrast.matrix=contrast.matrix)

#################################
# Part 4: Rerun for just predawn/noon
#################################
subdat<-which(info$Time %in% c("t0","t9"))
info<-info[subdat,]
counts<-counts[,subdat]
info$Time<-factor(info$Time, levels=c("t0", "t9"))
info$time_day_try<-factor(as.character(info$time_day_try), levels=c("w_d1_t0", "w_d1_t9", "d_d1_t0","d_d1_t9",
                                                                    "w_d2_t0", "w_d2_t9", "d_d2_t9","d_d2_t0"))
with(info, table(Time, Day, Treatment))

stats<-pipeLIMMA(counts=counts, info=info, block=info$id, formula="~ Treatment * Time + Day")
v<-stats$voom[["E"]]
stats.fullmodel.subtime<-stats$simpleStats
stats.allests.subtime<-stats$stats

f<-info$time_day_try
design <- model.matrix(~0+f)
contrast.matrix<-makeContrasts(fw_d1_t0-fd_d1_t0, fw_d1_t9-fd_d1_t9,
                               fw_d2_t0-fd_d2_t0, fw_d2_t9-fd_d2_t9,
                               levels=design)
lim.contrasts.subtime.all<-anovaLIMMA(counts=counts, design=design, block=info$id, contrast.matrix=contrast.matrix)

contrast.matrix<-makeContrasts(fw_d1_t0-fd_d1_t0,
                               fw_d2_t0-fd_d2_t0,
                               levels=design)
lim.contrasts.subtime.predawn<-anovaLIMMA(counts=counts, design=design, block=info$id, contrast.matrix=contrast.matrix)

contrast.matrix<-makeContrasts(fw_d1_t9-fd_d1_t9,
                               fw_d2_t9-fd_d2_t9,
                               levels=design)
lim.contrasts.subtime.midday<-anovaLIMMA(counts=counts, design=design, block=info$id, contrast.matrix=contrast.matrix)

save(stats.fullmodel, stats.allests, lim.contrasts, pca, v,
     stats.fullmodel.subtime, stats.allests.subtime,
     lim.contrasts.subtime.all, lim.contrasts.subtime.predawn, lim.contrasts.subtime.midday,
     file="/Users/John/Desktop/dropbox/Switchgrass_PlantPhys/stats_output/pots_allstats.RData")
