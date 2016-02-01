counts<-read.csv("./data/cylinder/goodCylinder_counts.csv")
info<-read.csv("./data/cylinder/goodCylinder_info.csv")
info$time<-factor(info$time, levels=c("dawn","midday", "dusk"))
info$trt<-factor(info$trt, levels=c("Wet","Dry"))

info$time_trt<-with(info, paste(time,trt, sep="_"))
info$time_trt<-factor(info$time_trt, levels=c("dawn_Dry", "dawn_Wet", "midday_Wet", "midday_Dry", "dusk_Dry","dusk_Wet"))
info$time_trt.mat<-as.numeric(info$time_trt)
save(info, counts, file="/Users/JLovell/Dropbox/Switchgrass_PlantPhys/physGenomicsPVFinal/data/cylinder.RData")
