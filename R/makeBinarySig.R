makeBinarySig<-function(x, alpha=0.05, what="q.value"){
  sig.q<-data.matrix(x[,grep(what, colnames(x))])
  if(length(grep(what, colnames(x)))==1){
    sig.q<-ifelse(sig.q<=alpha,1,0)
  }else{
    for(i in colnames(sig.q)) {
      sig.q[,i]<-ifelse(sig.q[,i]<=alpha,1,0)
      cat(i, sum(sig.q[,i]),"\n", sep="\t")
    }
  }
  return(sig.q)
}
