getEbayesStats<-function(fit,design){
  tests<-attr(design, "dimnames")[[2]][1:2]
  tests.out2<-lapply(tests, function(x){
    out2<-data.frame(fit$t[,grep(x, colnames(fit$t))],
                     fit$lods[,grep(x, colnames(fit$t))],
                     fit$p.value[,grep(x, colnames(fit$t))])
    colnames(out2)<-paste("ebayes",x,c("t","lods","p.value"),sep="_")
    n<-colnames(out2)[grep("p.value", names(out2))]
    out2[!is.finite(out2[,n]) | is.na(out2[,n]),n]<-1
    out2$qv<-qvalue(out2[,n])$qvalue
    colnames(out2)[which(colnames(out2)=="qv")]<-paste("ebayes",x,"q.value",sep="_")
    out2
  })
  tests.out2<-do.call(cbind, tests.out2)
  return(tests.out2)
}
