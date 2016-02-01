desStats<-function(...){
  res<-results(...)
  res$qvalue<-qvalue(res$pvalue)$qvalue
  par(mfrow=c(2,1))
  hist(res$pvalue, main=paste("pvalue",sum(res$pvalue<=0.05), sep="..."))
  hist(res$qvalue, main=paste("qvalue",sum(res$qvalue<=0.05), sep="..."))
  return(res)
}
