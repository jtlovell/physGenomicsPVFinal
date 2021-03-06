library(physGenomicsPVFinal)
rm(list=ls())
load("/Users/John/Desktop/dropbox/Switchgrass_PlantPhys/stats_output/cylinder_allstats.RData")
data(cylinder)

### for conviencence, rename info and counts
info<-infocyl
counts<-countscyl
library(ggplot2)
library(qdap)

# Summaries of genetic and physiological space

ggplot(info, aes(x=psi, y=RWC, col=time))+geom_point(aes(shape=trt), size=4) + theme_bw() +
  scale_shape_manual(values=c(2,19))+  scale_color_manual(values=c("cornflowerblue","darkred","green"))+
  scale_y_continuous("Leaf Relative Water Content (RWC)") + scale_x_continuous("Leaf Water Potential (MPa)") +
  ggtitle("Physiology data from Plants in cylinders")
ggplot(pca, aes(x=PC1, y=PC2, col=time))+geom_point(aes(shape=trt), size=4) + theme_bw()+
  scale_shape_manual(values=c(2,19))+  scale_color_manual(values=c("cornflowerblue","darkred","green"))+
  scale_y_continuous("PCA #2") + scale_x_continuous("PCA #1") +
  ggtitle("Gene Expression PCA")

#################################
# Part 1: Full model
# stats<-pipeLIMMA(counts=counts, info=info, block=info$pot, formula="~ trt * time")
# v<-stats$voom[["E"]]
# stats.fullmodel<-stats$simpleStats
# stats.allests<-stats$stats
#################################
pqHists(stats.allests, what.p="ebayes_trtdry_p.value",
        what.q="ebayes_trtdry_q.value",
        main="treatment effect", breaks=100)

ggplot(stats.allests, aes(x=ebayes_trtdry_coefficients, y=(-log10(ebayes_trtdry_p.value)), col=ebayes_trtdry_q.value<=0.05))+geom_point()+
  scale_color_manual(values=c(rgb(0,0,0,.5), rgb(1,0,0,.5)),guide = guide_legend(title = "is significant?"))+
  scale_y_continuous("-log10 P-value")+
  scale_x_continuous("Coefficient estimate of Treatment (Cylinder, all Times)")+
  theme_bw()+
  ggtitle("volcano plot for Cylinider data \n ~ Treatment * TimeOfDay")


#################################
# Part 2: All pairwise contrasts
# f<-info$time_trt
# design <- model.matrix(~0+f)
# contrast.matrix<-makeContrasts(fdawn_Wet-fdawn_Dry, fmidday_Wet-fmidday_Dry, fdusk_Wet-fdusk_Dry,
#                                levels=design)
# lim.contrasts<-anovaLIMMA(counts=counts, design=model.matrix(~0+f), block=info$pot, contrast.matrix=contrast.matrix)
#################################
pqHists(lim.contrasts, what.p="p.value_Ftest",
        what.q="q.value_Ftest",
        main="overall treatment effect", breaks=100)
pqHists(lim.contrasts, what.p="p.value_fmidday_Wet...fmidday_Dry",
        what.q="q.value_fmidday_Wet...fmidday_Dry",
        main="midday treatment effect", breaks=100)
pqHists(lim.contrasts, what.p="p.value_fdusk_Wet...fdusk_Dry",
        what.q="q.value_fdusk_Wet...fdusk_Dry",
        main="dusk treatment effect", breaks=100)
pqHists(lim.contrasts, what.p="q.value_fdawn_Wet...fdawn_Dry",
        what.q="q.value_fdawn_Wet...fdawn_Dry",
        main="dawn treatment effect", breaks=100)

sig.q.05<-makeBinarySig(lim.contrasts, what="q.value", alpha=0.05)
counts2Venn(x=sig.q.05, cols=c(2, 3, 4), names=c("dawn","midday","dusk"), colors=c("darkblue","cyan","red"),
            main="comparison of genes affected by \n ~Treatment at each time (alpha = 0.05)")
sig.q.1<-makeBinarySig(lim.contrasts, what="q.value", alpha=0.1)
counts2Venn(x=sig.q.1, cols=c(2, 3, 4), names=c("dawn","midday","dusk"), colors=c("darkblue","cyan","red"),
            main="comparison of genes affected by \n ~Treatment at each time (alpha = 0.1)")
sig.q.2<-makeBinarySig(lim.contrasts, what="q.value", alpha=0.2)
counts2Venn(x=sig.q.2, cols=c(2, 3, 4), names=c("dawn","midday","dusk"), colors=c("darkblue","cyan","red"),
            main="comparison of genes affected by \n ~Treatment at each time (alpha = 0.2)")

opar<-par()
info$time_trt_short<-multigsub(c("dawn","midday","dusk","Wet","Dry","_"),c("am","md","pm","w","d","."),as.character(info$time_trt))
v.means<-voom2MeanHeatMaps(v=v[stats.allests[,"ebayes_trtdry:timedusk_q.value"]<=0.05,], grps=info$time_trt_short, rowids=info$id,thresh=10)
par(opar)

ggplot(lim.contrasts, aes(x=coef_fdawn_Wet...fdawn_Dry, y=(-log10(p.value_fdawn_Wet...fdawn_Dry)), col=q.value_fdawn_Wet...fdawn_Dry<=0.05))+geom_point()+
  scale_color_manual(values=c(rgb(0,0,0,.5), rgb(1,0,0,.5)),guide = guide_legend(title = "is significant?"))+
  scale_y_continuous("-log10 P-value")+
  scale_x_continuous("Coefficient estimate of Treatment (Dawn sampling)")+
  theme_bw()+
  ggtitle("volcano plot for Cylinider data \n ~ Wet vs. Dry at Dawn")

ggplot(lim.contrasts, aes(x=coef_fmidday_Wet...fmidday_Dry, y=(-log10(p.value_fmidday_Wet...fmidday_Dry)), col=q.value_fmidday_Wet...fmidday_Dry<=0.05))+geom_point()+
  scale_color_manual(values=c(rgb(0,0,0,.5), rgb(1,0,0,.5)),guide = guide_legend(title = "is significant?"))+
  scale_y_continuous("-log10 P-value")+
  scale_x_continuous("Coefficient estimate of Treatment (midday sampling)")+
  theme_bw()+
  ggtitle("volcano plot for Cylinider data \n ~ Wet vs. Dry at midday")

ggplot(lim.contrasts, aes(x=coef_fdusk_Wet...fdusk_Dry, y=(-log10(p.value_fdusk_Wet...fdusk_Dry)), col=q.value_fdusk_Wet...fdusk_Dry<=0.05))+geom_point()+
  scale_color_manual(values=c(rgb(0,0,0,.5), rgb(1,0,0,.5)),guide = guide_legend(title = "is significant?"))+
  scale_y_continuous("-log10 P-value")+
  scale_x_continuous("Coefficient estimate of Treatment (dusk sampling)")+
  theme_bw()+
  ggtitle("volcano plot for Cylinider data \n ~ Wet vs. Dry at dusk")

stats.contrast<-lim.contrasts
stats.contrast$sig.duskdawn<-with(stats.contrast,ifelse(q.value_fdusk_Wet...fdusk_Dry<=0.05 & q.value_fdawn_Wet...fdawn_Dry<=0.05,"both",
                                                        ifelse(q.value_fdusk_Wet...fdusk_Dry<=0.05,"dusk",
                                                               ifelse(q.value_fdawn_Wet...fdawn_Dry<=0.05,"dawn","NS"))))
stats.contrast$sig.duskmidday<-with(stats.contrast,ifelse(q.value_fdusk_Wet...fdusk_Dry<=0.05 & q.value_fmidday_Wet...fmidday_Dry<=0.05,"both",
                                                          ifelse(q.value_fdusk_Wet...fdusk_Dry<=0.05,"dusk",
                                                                 ifelse(q.value_fmidday_Wet...fmidday_Dry<=0.05,"midday","NS"))))
stats.contrast$sig.middaydawn<-with(stats.contrast,ifelse(q.value_fmidday_Wet...fmidday_Dry<=0.05 & q.value_fdawn_Wet...fdawn_Dry<=0.05,"both",
                                                          ifelse(q.value_fmidday_Wet...fmidday_Dry<=0.05,"midday",
                                                                 ifelse(q.value_fdawn_Wet...fdawn_Dry<=0.05,"dawn","NS"))))

stats.contrast$sig.duskdawn<-factor(stats.contrast$sig.duskdawn, levels=c("both","dusk","dawn","NS"))
stats.contrast$sig.duskmidday<-factor(stats.contrast$sig.duskmidday, levels=c("both","dusk","midday","NS"))
stats.contrast$sig.middaydawn<-factor(stats.contrast$sig.middaydawn, levels=c("both","dawn","midday","NS"))

ggplot(stats.contrast, aes(x=t_fdusk_Wet...fdusk_Dry, y=t_fdawn_Wet...fdawn_Dry, col=sig.duskdawn))+geom_point()+
  scale_color_manual(values=c("purple", "red", "skyblue", rgb(0,0,0,.5)),guide = guide_legend(title = "is significant?"))+
  scale_y_continuous("T statistic of Treatment (dawn sampling)")+
  scale_x_continuous("T statistic of Treatment (dusk sampling)")+
  theme_bw()+
  ggtitle("Time correlations for Cylinider data \n ~ Wet vs. Dry, Dusk and Dawn")

ggplot(stats.contrast, aes(x=t_fdusk_Wet...fdusk_Dry, y=t_fmidday_Wet...fmidday_Dry, col=sig.duskmidday))+geom_point()+
  scale_color_manual(values=c("purple", "red", "skyblue", rgb(0,0,0,.5)),guide = guide_legend(title = "is significant?"))+
  scale_y_continuous("T statistic of Treatment (midday sampling)")+
  scale_x_continuous("T statistic of Treatment (dusk sampling)")+
  theme_bw()+
  ggtitle("Time correlations for Cylinider data \n ~ Wet vs. Dry, Midday and Dusk")

ggplot(stats.contrast, aes(x=t_fmidday_Wet...fmidday_Dry, y=t_fdawn_Wet...fdawn_Dry, col=sig.middaydawn))+geom_point()+
  scale_color_manual(values=c("purple", "red", "skyblue", rgb(0,0,0,.5)),guide = guide_legend(title = "is significant?"))+
  scale_y_continuous("T statistic of Treatment (dawn sampling)")+
  scale_x_continuous("T statistic of Treatment (midday sampling)")+
  theme_bw()+
  ggtitle("Time correlations for Cylinider data \n ~ Wet vs. Dry, Midday and Dawn")
