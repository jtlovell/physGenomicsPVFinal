hhmmss2Time<-function(x){
  x<-as.character(x)
  h<-as.numeric(sapply(x, function(y) strsplit(y, ":")[[1]][1]))
  m<-as.numeric(sapply(x, function(y) strsplit(y, ":")[[1]][2]))
  s<-as.numeric(sapply(x, function(y) strsplit(y, ":")[[1]][3]))
  d<-data.frame(h,m,s)
  d$time<-with(d, (h*60)+m+(s/60))
  d$time<-round(d$time-(min(d$time)),2)
  return(d$time)
}
