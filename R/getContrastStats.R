getContrastStats<-function(fit,contrasts, names){
  pvals<-fit$p.value
  ts<-fit$t
  lods<-fit$lods
  tests.out2<-lapply(1:length(names), function(x){
    out2<-data.frame(fit$t[,x],
                     fit$lods[,x],
                     fit$p.value[,x])
    colnames(out2)<-paste("ebayes",names[x],c("t","lods","p.value"),sep="_")
    n<-colnames(out2)[grep("p.value", names(out2))]
    out2[!is.finite(out2[,n]) | is.na(out2[,n]),n]<-1
    out2$qv<-qvalue(out2[,n])$qvalue
    colnames(out2)[which(colnames(out2)=="qv")]<-paste("ebayes",names[x],"q.value",sep="_")
    out2
  })
  tests.out2<-do.call(cbind, tests.out2)
  return(tests.out2)
}
