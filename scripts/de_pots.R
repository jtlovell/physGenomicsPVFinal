library(devtools)
install_github("jtlovell/physGenomicsPVFinal")
library(physGenomicsPVFinal)
rm(list=ls())
data(potsGreenhouse)

info<-info.wetdry
counts<-counts.recovery

#################################
# Part 4.1: Full model
#################################
stats<-pipeLIMMA(counts=counts.wetdry, info=info.wetdry, block=info$id, formula="~ Treatment * Time + Day")
voom.trtByMonth<-stats$voom[["E"]]
stats.fullmodel<-stats$simpleStats
stats.allests<-stats$stats

#################################
# Part 4.2: Run PCA
#################################
pca<-voom2PCA(v=voom.trtByMonth, info=info, ids=info$id)
ggplot(pca, aes(x=PC1, y=PC2, col=Day, shape=Treatment))+theme_bw() +facet_grid(Time~Treatment)+geom_point(size=4)+scale_shape_manual(values=c(2,19))+
  geom_text(aes(label=Sample))


