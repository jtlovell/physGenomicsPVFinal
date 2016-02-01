limmaStats<-function(x, what="q.value", threshold=0.05){
  out<-data.frame(gene=x[,1])
  for(i in names(x)[grep(what,names(x))]){
    j=i
    cat(i, sum(x[,i]<=threshold), "\n")
    temp<-data.frame(x[,i])
    names(temp)<-i
    #hist(temp[,i], breaks=200, main=j, xlab=j)
    out<-data.frame(out,temp)
  }
  return(out)
}

