library(physGenomicsPVFinal)
rm(list=ls())

load("/Users/John/Desktop/dropbox/Switchgrass_PlantPhys/stats_output/shelter201314_allstats.RData")
data(tmpwfc20132014_2treatments)

### for conviencence, rename info and counts
info<-info1314
counts<-counts1314
library(ggplot2)
#################################
### Part 1: Statistical presentation of main model
# stats<-pipeLIMMA(counts=counts, info=info, block=info$sb_unique, formula="~ Treatment * Location + Year + order")
# v<-stats$voom[["E"]]
# stats.fullmodel<-stats$simpleStats
# stats.allests<-stats$stats
#################################
ggplot(info, aes(x=PDWP, y=MDWP, col=Location, shape=Treatment))+geom_point(size=4)+ theme_bw() +
  stat_smooth(method="lm", se=F, lty=2, alpha=.5, lwd=.5)+
  scale_color_manual(values=c("darkblue","orange"))+
  scale_shape_manual(values=c(2,19))+
  facet_wrap(~Year)+
  scale_x_continuous("Pre-dawn Leaf Water Potential (MPa)")+
  scale_y_continuous("Mid Day Leaf Water Potential (MPa)")+
  ggtitle("Effect of sampling order on measured leaf water potential")


pqHists(stats.fullmodel, what.p="pvalue", what.q="qvalue", main="full model", breaks=100)
pqHists(stats.allests, what.p="ebayes_treatmentlow:locationwfc_p.value",
        what.q="ebayes_treatmentlow:locationwfc_q.value",
        main="location*treatment effect", breaks=100)
pqHists(stats.allests, what.p="ebayes_order_p.value",
        what.q="ebayes_order_q.value",
        main="sampling order effect", breaks=100)
pqHists(stats.allests, what.p="ebayes_year2014_p.value",
        what.q="ebayes_year2014_q.value",
        main="year effect", breaks=100)
pqHists(stats.allests, what.p="ebayes_locationwfc_p.value",
        what.q="ebayes_locationwfc_q.value",
        main="location effect", breaks=100)
pqHists(stats.allests, what.p="ebayes_treatmentlow_p.value",
        what.q="ebayes_treatmentlow_q.value",
        main="treatment effect", breaks=100)

sig.q.05<-makeBinarySig(stats.allests, what="q.value", alpha=0.05)
counts2Venn(x=sig.q.05, cols=c(2, 3, 4), names=c("trt","loc","year"), colors=c("darkblue","cyan","red"),
            main="comparison of genes affected by \n ~Treatment * Location + Year  in 2013/2014")
counts2Venn(x=sig.q.05, cols=c(2, 3, 6), names=c("trt","loc","trt*location"), colors=c("darkblue","cyan","red"),
            main="comparison of genes affected by \n ~Treatment * Location  in 2013/2014")

opar<-par()
v.means<-voom2MeanHeatMaps(v=v[stats.allests[,"ebayes_treatmentlow:locationwfc_p.value"]<=0.001,], grps=info$Year_Location_Treatment,rowids=info$ID,thresh=8)
par(opar)
#################################
### Part 2: Comparison of contrasts
# design <- model.matrix(~ 0+factor(Year_Location_Treatment.mat) + order , data=info)
# colnames(design)[1:8]<-c("TMP_13_high","TMP_13_low","WFC_13_high","WFC_13_low",
#                          "TMP_14_high","TMP_14_low","WFC_14_high","WFC_14_low")
# contrast.matrix <- makeContrasts(TMP_13_high-TMP_13_low, WFC_13_high-WFC_13_low, TMP_14_high-TMP_14_low, WFC_14_high-WFC_14_low,
#                                  levels=design)
# lim.contrasts<-anovaLIMMA(counts=counts, design=design, block=info$Sub_Block, contrast.matrix=contrast.matrix)
#################################
sig.q.05<-makeBinarySig(lim.contrasts, what="q.value", alpha=0.05)
counts2Venn(x=sig.q.05, cols=c(2, 4, 5), names=c("tmp13","tmp14","wfc14"), colors=c("darkblue","cyan","red"),
            main="comparison of genes affected by \n ~ Treatment Contrasts within Year*Site \n alpha=0.05")

sig.q.1<-makeBinarySig(lim.contrasts, what="q.value", alpha=0.1)
counts2Venn(x=sig.q.1, cols=c(2, 4, 5), names=c("tmp13","tmp14","wfc14"), colors=c("darkblue","cyan","red"),
            main="comparison of genes affected by \n ~ Treatment Contrasts within Year*Site \n alpha=0.1")

sig.q.2<-makeBinarySig(lim.contrasts, what="q.value", alpha=0.2)
counts2Venn(x=sig.q.2, cols=c(2, 4, 5), names=c("tmp13","tmp14","wfc14"), colors=c("darkblue","cyan","red"),
            main="comparison of genes affected by \n ~ Treatment Contrasts within Year*Site \n alpha=0.2")

stats.cnt<-data.frame(lim.contrasts)
stats.cnt$sig<-with(stats.cnt, ifelse(q.value_TMP_14_high...TMP_14_low<=0.05 & q.value_TMP_13_high...TMP_13_low <= 0.05,"sigBoth",
                                             ifelse(q.value_TMP_14_high...TMP_14_low<=0.05,"sig2013",
                                                    ifelse(q.value_TMP_13_high...TMP_13_low<=0.05,"sig2014","NS"))))
stats.cnt$sig<-factor(stats.cnt$sig, levels=c("sigBoth","sig2013","sig2014","NS"))

ggplot(stats.cnt, aes(x=t_TMP_13_high...TMP_13_low, y=t_TMP_14_high...TMP_14_low, color=sig))+  theme_bw()+
  geom_point()+scale_color_manual(values=c(rgb(0,1,0,.5), rgb(1,0,0,.5), rgb(0,0,1,.5), rgb(0,0,0,.2)),guide = guide_legend(title = "is significant?"))+
  theme_bw()+geom_hline(yintercept=0, lty=2, lwd=.5, col="grey")+ geom_vline(xintercept=0, lty=2, lwd=.5, col="grey")+
  scale_y_continuous("t estimate (Wet vs. Dry) in Temple, 2014", limits=c(min(c(stats.cnt$t_TMP_13_high...TMP_13_low,stats.cnt$t_TMP_14_high...TMP_14_low)),max(c(stats.cnt$t_TMP_13_high...TMP_13_low,stats.cnt$t_TMP_14_high...TMP_14_low))))+
  scale_x_continuous("t estimate (Wet vs. Dry) in Temple, 2013", limits=c(min(c(stats.cnt$t_TMP_13_high...TMP_13_low,stats.cnt$t_TMP_14_high...TMP_14_low)),max(c(stats.cnt$t_TMP_13_high...TMP_13_low,stats.cnt$t_TMP_14_high...TMP_14_low))))+
  theme(panel.grid.major = element_blank() ,  panel.grid.minor = element_blank())+
  ggtitle("comparison of effects in Temple-specific treatment contrasts")

stats.cnt$id<-rownames(stats.cnt)
stats.cnt.tmp<-melt(stats.cnt, id.vars=c("id","sig"), measure.var=c("t_TMP_14_high...TMP_14_low", "t_TMP_13_high...TMP_13_low"))
stats.cnt.tmp$variable<-factor(stats.cnt.tmp$variable, levels=c("t_TMP_13_high...TMP_13_low","t_TMP_14_high...TMP_14_low"))

ggplot(stats.cnt.tmp[stats.cnt$sig!="NS",], aes(x=variable, y=value, group=id, color=sig))+  theme_bw()+
  geom_line()+scale_color_manual(values=c(rgb(0,1,0,.5), rgb(1,0,0,.5), rgb(0,0,1,.5)),guide = guide_legend(title = "is significant?"))+
  scale_x_discrete("Experimental Year",labels=c("2013","2014"),expand=c(.1,.1))+
  scale_y_continuous("t-statistic estimate (Wet vs. Dry) in Temple")+
  theme(panel.grid.major = element_blank() ,  panel.grid.minor = element_blank())+
  ggtitle("comparison of effects in Temple-specific treatment contrasts \n only significantly DE genes")


#################################
# Part 3: Plotting of PCAs
#################################
library(ggplot2)

ggplot(pca, aes(x=PC1, y=PC2, col=Location, shape=Treatment))+geom_point(size=4)+ theme_bw() +
  scale_color_manual(values=c("darkblue","orange"))+
  scale_shape_manual(values=c(2,19))+
  facet_wrap(~Year)+
  scale_y_continuous("PCA #2") + scale_x_continuous("PCA #1") +  ggtitle("Normalized expresssion PCA")

ggplot(pca, aes(x=MDWP, y=PC1, col=Location, shape=Treatment))+geom_point(size=4)+ theme_bw() +
  scale_color_manual(values=c("darkblue","orange"))+
  scale_shape_manual(values=c(2,19))+
  facet_wrap(~Year, scales="free")+
  scale_y_continuous("PCA #1") + scale_x_continuous("Mid Day Leaf Water Potential (MPa)") +  ggtitle("Normalized expresssion PCA")
ggplot(pca, aes(x=MDWP, y=PC2, col=Location, shape=Treatment))+geom_point(size=4)+ theme_bw() +
  scale_color_manual(values=c("darkblue","orange"))+
  scale_shape_manual(values=c(2,19))+
  facet_wrap(~Year, scales="free")+
  scale_y_continuous("PCA #2") + scale_x_continuous("Mid Day Leaf Water Potential (MPa)") +  ggtitle("Normalized expresssion PCA")

#################################
# Part 4: Effect of midday water potential
# stats<-pipeLIMMA(counts=counts, info=info, block=info$sb_unique, formula="~ MDWP * Location + Year + order")
# stats.fullmodel.mdwp<-stats$simpleStats
# stats.allests.mdwp<-stats$stats
# ## Run model with MDWP as the predictor, without controlling for location or year
# stats<-pipeLIMMA(counts=counts, info=info, block=info$sb_unique, formula="~ MDWP + order")
# stats.fullmodel.mdwponly<-stats$simpleStats
# stats.allests.mdwponly<-stats$stats
#################################
mdwp.q.05<-makeBinarySig(stats.allests.mdwponly, what="mdwp_q.value", alpha=0.05)
trt.q.05<-makeBinarySig(lim.contrasts, what="q.value_Ftest", alpha=0.05)
mdwp.q.1<-makeBinarySig(stats.allests.mdwponly, what="mdwp_q.value", alpha=0.1)
trt.q.1<-makeBinarySig(lim.contrasts, what="q.value_Ftest", alpha=0.1)
trt.mdwp.qs<-data.frame(trt.05=trt.q.05,mdwp.05=mdwp.q.05, trt.1=trt.q.1, mdwp.1=mdwp.q.1)
colSums(trt.mdwp.qs)
par(mfrow=c(1,1))
counts2Venn(x=trt.mdwp.qs, cols=c(1,2), names=c("trt","MDWP"), colors=c("darkblue","red"),
            main="comparison of  genes affected by \n ~ Treatment vs. Midday Water Potential F-tests alpha = 0.05")
counts2Venn(x=trt.mdwp.qs, cols=c(3,4), names=c("trt","MDWP"), colors=c("darkblue","red"),
            main="comparison of  genes affected by \n ~ Treatment vs. Midday Water Potential F-tests alpha = 0.1")
