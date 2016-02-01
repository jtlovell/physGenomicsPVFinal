pipeLIMMA<-function(counts, info, formula, use.qualityWeights=TRUE, block, tests="all",geneIDs=NA, useBlock=TRUE, getTopTable=TRUE, getEbayes=T, contrasts=NULL, ...){
  if(is.na(geneIDs)){
    geneIDs<-rownames(counts)
  }

  design<-model.matrix(as.formula(formula), data = info)
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

  fit<-eBayes(fit)
  out<-data.frame(gene=geneIDs,
                  sigma=fit$sigma,
                  s2.post=fit$s2.post,
                  Amean=fit$Amean)
  if(tests=="all"){
    tests<-attr(design, "dimnames")[[2]]
  }
  tests.out<-lapply(tests, function(x){
    if(getEbayes & getTopTable){
      out2<-data.frame(fit$stdev.unscaled[,x],
                       fit$coefficients[,x],
                       fit$lods[,x],
                       fit$p.value[,x],
                       qvalue(fit$p.value[,x])$qvalue)
      colnames(out2)<-paste("ebayes",x,c("stdev.unscaled","coefficients","lods","p.value","q.value"),sep="_")
      out3<-data.frame(toptable(fit, p.value=1, coef=x, number=100000))
      out3<-out3[,c("logFC","t","B")]
      colnames(out3)<-paste("tt",x,colnames(out3),sep="_")
      out2<-data.frame(out2, out3)
    }else{
      if(getTopTable){
        out2<-data.frame(toptable(fit, p.value=1, coef=x, number=100000))
        out3<-out3[,c("logFC","t","B")]
        colnames(out2)<-paste("tt", x,colnames(out2),sep="_")
      }else{
        out2<-data.frame(fit$stdev.unscaled[,x],
                        fit$coefficients[,x],
                        fit$lods[,x],
                        fit$p.value[,x],
                        qvalue(fit$p.value[,x])$qvalue)
        colnames(out2)<-paste("ebayes",x,c("stdev.unscaled","coefficients","lods","p.value","q.value"),sep="_")
      }
    }
    out2
  })
  tests.out2<-do.call(cbind, tests.out)
  all.out<-cbind(data.frame(out),tests.out2)
  colnames(all.out)<-tolower(colnames(all.out))
  return(list(stats=all.out, voom=v, lmfit=fit, countsSize=y))
}
