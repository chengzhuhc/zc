c2KEGG <- getGmt("markes.gmt",geneIdType=SymbolIdentifier())
ssGSEA.result <- gsva(as.matrix(data.cut), c2KEGG,method='ssgsea', min.sz=3, max.sz=500, verbose=TRUE)
write.csv(ssGSEA.result, paste0(analysis.dir,"../markes_output.csv"), quote = F)

cor_point=function(x,y,method='Pearson',top_col='#D55E00',right_col='#009E73'
                   ,ylab='y expression',xlab='x expression',title=NULL
                   ,marginal.type=c("histogram", "boxplot", "density", "violin", "densigram")[1]){
  library(ggstatsplot)
  dat=data.frame(X=x,Y=y)
  tp='nonparametric'
  if(method=='Pearson'|method=='pearson'){
    tp='parametric'
  }
  g1=ggscatterstats(data = dat, 
                    x = X, 
                    y = Y
                    ,type = tp
                    ,xfill = top_col
                    ,yfill = right_col
                    ,xlab = xlab
                    ,ylab=ylab
                    ,marginal.type = marginal.type
                    ,title = title)
  return(g1)  
}

plotcor=list()
for (aaa in 1:nrow(ssGSEA.result)) {
    aaa1=rownames(ssGSEA.result)[aaa]
    plotcor[[aaa1]]=cor_point(x=as.numeric(data.cut[genes,]),y=as.numeric(ssGSEA.result[aaa1,]),top_col=mypal[1],right_col=mypal[2]
                            ,ylab=aaa1,xlab=paste0('Log2 (',genes,' TPM + 1)'),marginal.type='density',method = 'nonparametric')
      
}

