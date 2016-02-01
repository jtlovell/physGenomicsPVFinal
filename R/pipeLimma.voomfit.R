pipeLIMMA.voomfit<-function(counts, info, design, use.qualityWeights=TRUE, block, tests="all",geneIDs=NA, useBlock=TRUE, getTopTable=TRUE, getEbayes=T, contrasts=NULL, ...){
  if(is.na(geneIDs)){
    geneIDs<-rownames(counts)
  }
  y <- DGEList(counts = counts)
  y <- calcNormFactors(y)

  if(use.qualityWeights){
    v <- voomWithQualityWeights(y, design=design, plot = F)
  }else{
    v <- voom(y, design=design, plot = F)
  }
  if(useBlock){
    dupcor <- duplicateCorrelation(counts,design, block=as.factor(block))
    fit <- lmFit(v, design=design, correlation=dupcor$consensus, block=as.factor(block))
  }else{
    fit <- lmFit(v, design=design)
  }
  efit<-eBayes(fit)
  return(list(normFactors=y, voom=v, lmFit=fit, ebayes=efit))
}
