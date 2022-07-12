load('TCGA_XXX.RData')
tcga.data=log(dt+1)
load('TCGA_XXX.RData')
gtex.data=log(dt+1)
load('Clinical.RData')
gene.type=read.csv('GeneTag.txt',stringsAsFactors = F,row.names = 1,sep = '\t',check.names = F)
gene.type1=gene.type[which(gene.type$TYPE=='protein_coding'),]
nor.data=tcga.data[,grep('-11',colnames(tcga.data))]
comsamples=intersect(c(paste0(rownames(clinical),'-01'),paste0(rownames(clinical),'-03'),paste0(rownames(clinical),'-06')),colnames(tcga.data))
tum.data=tcga.data[,comsamples]
clin.cut=clinical[substr(comsamples,1,12),]
comgenes=intersect(intersect(rownames(tcga.data),rownames(gtex.data)),gene.type1$SYMBOL)
datas=cbind.data.frame(tum.data[comgenes,],nor.data[comgenes,],gtex.data[comgenes,])
gs.file=c(rep('Tumor',ncol(tum.data)),rep('Normal',ncol(nor.data)),rep('Normal',ncol(gtex.data)))
datas=datas[which(apply(datas,1,function(x){return(sum(x>0))})>0.5*ncol(datas)),]

limma_DEG=function(exp,group,ulab,dlab){
  library(limma)
  ind1=which(group==ulab)
  ind2=which(group==dlab)
  
  sml <- c(rep('G1',length(ind1)),rep('G0',length(ind2)))    # set group names
  eset=exp[,c(ind1,ind2)]
  fl <- as.factor(sml)
  
  design <- model.matrix(~fl+0)
  colnames(design) <- levels(fl)
  cont.matrix<-makeContrasts(contrasts='G1-G0',levels=design)
  #print(head(eset))
  fit<-lmFit (eset,design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  #print(sml)
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(eset))
  return(list(Exp=eset,Group=group[c(ind1,ind2)],DEG=tT))
}

deg=limma_DEG(datas,gs.file,'Tumor','Normal')    
deg.result=deg$DEG[which(abs(as.numeric(deg$DEG$logFC)) > 1&as.numeric(deg$DEG$FDR) < 0.05),]

up.genes=rownames(deg.result)[which(deg.result$logFC > 0)]
down.genes=rownames(deg.result)[which(deg.result$logFC < 0)]

Volcano=function(logfc,pvalue,symbol=NULL,cutFC=1,cutPvalue=0.05
                    ,showText=NULL
                    ,colors=c(mypal[2],'grey',mypal[1])
                    ,xlim=NULL,ylim=NULL
                    ,legend.pos='tl'
                    ,ylab='-log10(FDR)',leg='',xlab='log2(FoldChange)'){
  library(ggplot2)
  pos=c(0,0)
  if(is.null(legend.pos)){
    pos='none'
  }else if(legend.pos=='tr'){
    pos=c(1,1)
  }else if(legend.pos=='br'){
    pos=c(1,0)
  }else if(legend.pos=='tl'){
    pos=c(0,1)
  }else if(legend.pos=='bl'){
    pos=c(0,0)
  }else{
    pos='right'
  }
  cange=rep('None',length(logfc))
  cange[which(logfc>cutFC&pvalue<cutPvalue)]='Up-regulation'
  cange[which(logfc< -cutFC&pvalue<cutPvalue)]='Down-regulation'
  if(is.null(symbol)){
    symbol=rep('',length(logfc))
    showText=NULL
  }
  vo.input=data.frame(logFC=logfc,FDR=pvalue,change=cange,SYMBOL=symbol)
  #print(head(vo.input))
  p1 <- ggplot(data = vo.input, 
               aes(x = logFC, 
                   y = -log10(FDR)))
  if (ylab == 'FDR') {
    p1=p1+geom_point(alpha=0.4, size=3.5, aes(color=change))+labs(x=bquote(~Log[2]~"(fold change)"), y=bquote(~-Log[10]~italic('FDR')), title="")
  }else{
    p1=p1+geom_point(alpha=0.4, size=3.5, aes(color=change))+labs(x=bquote(~Log[2]~"(fold change)"), y=bquote(~-Log[10]~italic('P-value')), title="")
  }
  p1=p1+scale_color_manual(values=colors,limits = c("Down-regulation",'None', "Up-regulation"),name=leg)+ scale_x_continuous(
    breaks = c(-10, -5, -cutFC, 0, cutFC, 5, 10), #刻度线的位置
    labels = c(-10, -5, -cutFC, 0, cutFC, 5, 10)) #x轴范围，两侧对称才好看
  p1=p1+geom_vline(xintercept=c(-cutFC,cutFC),lty=4,col="black",lwd=0.8)  
  p1=p1+geom_hline(yintercept = -log10(cutPvalue),lty=4,col="black",lwd=0.8)  
  p1=p1+theme_bw()
  p1=p1+theme(
    axis.text.y=element_text(family="serif",face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
    axis.title.y=element_text(family="serif",face="plain"), #设置y轴标题的字体属性
    legend.text=element_text(family="serif",face="plain", colour="black"  #设置图例的子标题的字体属性
    ),
    legend.title=element_text(family="serif",face="plain", colour="black" #设置图例的总标题的字体属性
    ),
    legend.justification=pos, legend.position=pos
    ,legend.background = element_rect(fill = NA, colour = NA)
  )
  if(is.null(showText)|is.null(symbol)){
    showText=c()
  }
  
  if(length(showText)>0){
    for_label <-vo.input[match(intersect(showText,vo.input$SYMBOL),vo.input$SYMBOL),]
    p1=p1+geom_point(size = 3, shape = 1, data = for_label)+ggrepel::geom_label_repel(
      aes(label = SYMBOL),
      data = for_label,
      color="black"
    )
  }
  if(!is.null(ylim)){
    p1=p1+ylim(ylim)
  }
  if(!is.null(xlim)){
    p1=p1+xlim(xlim)
  }
  p1=p1+theme(text=element_text(size=12,family="serif"))
  return(p1)
}


gene_volcano=Volcano(as.numeric(deg$DEG$logFC),as.numeric(deg$DEG$FDR),cutFC = cut.fc,colors=c('red','grey','blue')
                        ,cutPvalue = cut.p,symbol = rownames(deg$DEG),ylab = 'P-value',showText = symlist)

if (length(up.genes)>100){
  genes1=rownames(deg.dif.result.final)[order(deg.dif.result.final$logFC,decreasing = T)][1:100]
}else{
  genes1=up.genes
}
if (length(down.genes)>100){
  genes2=rownames(deg.dif.result.final)[order(deg.dif.result.final$logFC)][1:100]
}else{
  genes2=down.genes
}

indata=indata[c(genes1,genes2),]
indata=as.data.frame(remove_NA(indata))

annCol=as.data.frame(group)
rownames(annCol)=colnames(datas)
###################################
annColors <- list()
if (length(clinType)>1) {
  annColors[["group"]] <- c("G1"="#EA6767","G2"="#70C17A")
}else{
  annColors[["group"]] <- c("G1"="#EA6767","Normal"="#70C17A")
}
blank <- "   "
add.label <- str_pad(rep("   ",nrow(indata)), # 固定行名宽度并再右侧补齐" "
                     max(nchar(paste0(rownames(indata)))), 
                     side = "right")

gene_map=pheatmap(mat = indata
                  ,clustering_distance_rows = 'correlation'
                  ,clustering_distance_cols = 'correlation'
                  ,scale = "none" # 不标准化
                  ,border_color = NA
                  ,color = colorRampPalette(c(mypal1[2], "white", mypal1[1]))(100) # 例文配色
                  ,cluster_cols = T # 列不聚类
                  ,cluster_rows = T # 行不聚类
                  ,show_rownames = T # 不显示行名
                  ,show_colnames = F # 不显示列名
                  ,annotation_col = annCol # 列注释（注意要根据TIM-3的表达顺序排列
                  #,annotation_row = annrow
                  ,annotation_colors = annColors
                  ,labels_row = paste(add.label,sep=blank))

fun_clusterProfiler=function(genes,minGSSize=10,maxGSSize = 500, pAdjustMethod = "BH",pvalueCutoff=0.05,fdrCutoff = 0.5,qvalueCutoff = 0.5){
  library(org.Hs.eg.db)
  #pathways=readMatrix('/opt/shengxin/app/pathway_gids.txt',header=F,row=F)
  pathways=readMatrix('./pathway_gids.txt',header=F,row=F)
  eid=AnnotationDbi::select(org.Hs.eg.db,keys = as.character(pathways[,2]),keytype = 'ENTREZID',columns = c('SYMBOL'))
  pathways[,2]=eid$SYMBOL[match(pathways[,2],eid$ENTREZID)]
  pathways=pathways[!is.na(pathways[,2]),]
  
  #hg.tab=select(org.Hs.eg.db, keys=keys(org.Hs.eg.db,keytype = "ENTREZID"),columns=c('GO','GOALL')
  #       , keytype="ENTREZID")
  pCut=pvalueCutoff
  pvalueCutoff = fdrCutoff
  #qvalueCutoff = 0.5
  #print('Starting KEGG')
  kegg=clusterProfiler::enricher(genes, pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, 
                                 minGSSize = minGSSize, maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff
                                 ,TERM2GENE=data.frame(term =pathways[,1],gene=pathways[,2]),
                                 TERM2NAME = data.frame(term =pathways[,1],name=pathways[,3]))
  #print('Succ KEGG,Start GO_MF')
  eid=AnnotationDbi::select(org.Hs.eg.db,keys = as.character(genes),keytype = 'SYMBOL',columns = c('ENTREZID'))
  genes1=eid$ENTREZID
  genes1=as.character(genes1[!is.na(genes)])
  
  go.mf=clusterProfiler::enrichGO(genes1, org.Hs.eg.db, keyType = "ENTREZID", ont = "MF",
                                  pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, 
                                  qvalueCutoff = qvalueCutoff, minGSSize = minGSSize, maxGSSize = maxGSSize,
                                  readable = T, pool = FALSE)
  #print('Succ GO_MF,Start GO_CC')
  
  go.cc=clusterProfiler::enrichGO(genes1, org.Hs.eg.db, keyType = "ENTREZID", ont = "CC",
                                  pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, 
                                  qvalueCutoff = qvalueCutoff, minGSSize = minGSSize, maxGSSize = maxGSSize,
                                  readable = T, pool = FALSE)
  #print('Succ GO_CC,Start GO_BP')
  go.bp=clusterProfiler::enrichGO(genes1, org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                                  pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, 
                                  qvalueCutoff = qvalueCutoff, minGSSize = minGSSize, maxGSSize = maxGSSize,
                                  readable = T, pool = FALSE)
  enrich_tab=rbind(cbind((kegg@result),DB=rep('pathway_KEGG',nrow(kegg@result)))
                   ,cbind((go.bp@result),DB=rep('geneontology_Biological_Process',nrow(go.bp@result)))
                   ,cbind((go.cc@result),DB=rep('geneontology_Cellular_Component',nrow(go.cc@result)))
                   ,cbind((go.mf@result),DB=rep('geneontology_Molecular_Function',nrow(go.mf@result)))
  )
  enrich_tab=crbind2DataFrame(enrich_tab)
  colnames(enrich_tab)[c(2,5,6,9)]=c('description','pValue','FDR','size')
  enrich_tab$enrichmentRatio=unlist(lapply(strsplit(enrich_tab[,3],'/'), function(x){
    x1=as.numeric(x)
    return(x1[1]/x1[2])
  }))
  enrich_tab=enrich_tab[which(enrich_tab[,5]<pCut&enrich_tab[,6]<fdrCutoff&enrich_tab[,7]<pvalueCutoff),]
  return(list(KEGG=kegg,GO_BP=go.bp,GO_CC=go.cc,GO_MF=go.mf,Enrich_tab=enrich_tab))
}
dotplot_batch=function(pathway_data,db,colors,top=20,FDR=T){
  library(ggplot2)
  #library(ggsci)
  #db='pathway_KEGG'
  #pathway_data=up.result$KEGG@result
  if(nrow(pathway_data)>0){
    if(nrow(pathway_data)>top){
      pathway_data=pathway_data[order(pathway_data$p.adjust)[1:top],]
      tl=paste0(db)
    }else{
      tl=paste0(db)
    }
    desc=pathway_data$Description
    ndesc=c()
    for(de in desc){
      if(nchar(de)>50&length(grep(' ',de))>0){
        de1=unlist(strsplit(de,' '))
        d2=paste0(de1[(ceiling(length(de1)/2)+1):length(de1)],collapse = ' ')
        if(nchar(d2)>50){
          d2=paste0(substr(d2,0,47),'...')
        }
        de2=paste0(paste0(de1[1:ceiling(length(de1)/2)],collapse = ' '),'\n'
                   ,d2)
        ndesc=c(ndesc,de2)
      }else{
        ndesc=c(ndesc,de)
      }
    }
    pathway_data$Description=ndesc      
    pathway_data$pvalue[pathway_data$pvalue<1e-16]=1e-16
    pathway_data$p.adjust[pathway_data$p.adjust<1e-16]=1e-16
    
    for (i in 1:nrow(pathway_data)) {
      can11=as.data.frame(strsplit(pathway_data$GeneRatio[i],"/"),stringsAsFactors = FALSE)
      pathway_data$size[i]=round(as.numeric(can11[1,])/as.numeric(can11[2,]),2)
    }
    pathway_data=pathway_data[order(pathway_data$size),]
    
    bubble=ggplot(data = pathway_data, aes(x = size, y = Description)) 
    
    bubble=bubble+xlab('Enrichment Ratio')
    if(FDR){
      bubble1=bubble+geom_point(aes(size = Count,color = -log10(p.adjust))) + scale_color_gradient(low = ggsci::pal_npg()(2)[2], high = ggsci::pal_npg()(2)[1])
      bubble2=bubble+geom_point(aes(size = Count,color = -log10(p.adjust))) + scale_color_gradient(low = ggsci::pal_nejm()(2)[2], high = ggsci::pal_nejm()(2)[1])
      bubble3=bubble+geom_point(aes(size = Count,color = -log10(p.adjust))) + scale_color_gradient(low = ggsci::pal_lancet()(2)[1], high = ggsci::pal_lancet()(2)[2])
      bubble4=bubble+geom_point(aes(size = Count,color = -log10(p.adjust))) + scale_color_gradient(low = ggsci::pal_jama()(2)[1], high = ggsci::pal_jama()(2)[2])
      bubble5=bubble+geom_point(aes(size = Count,color = -log10(p.adjust))) + scale_color_gradient(low = ggsci::pal_jco()(2)[1], high = ggsci::pal_jco()(2)[2])
    }else{
      bubble1=bubble+geom_point(aes(size = Count,color = -log10(pValue))) + scale_color_gradient(low = ggsci::pal_npg()(2)[2], high = ggsci::pal_npg()(2)[1])
      bubble2=bubble+geom_point(aes(size = Count,color = -log10(pValue))) + scale_color_gradient(low = ggsci::pal_nejm()(2)[2], high = ggsci::pal_nejm()(2)[1])
      bubble3=bubble+geom_point(aes(size = Count,color = -log10(pValue))) + scale_color_gradient(low = ggsci::pal_lancet()(2)[1], high = ggsci::pal_lancet()(2)[2])
      bubble4=bubble+geom_point(aes(size = Count,color = -log10(pValue))) + scale_color_gradient(low = ggsci::pal_jama()(2)[1], high = ggsci::pal_jama()(2)[2])
      bubble5=bubble+geom_point(aes(size = Count,color = -log10(pValue))) + scale_color_gradient(low = ggsci::pal_jco()(2)[1], high = ggsci::pal_jco()(2)[2])
    }
    if (colors=='npg') {
      bubble=bubble1+ggsci::scale_fill_npg()+ggplot2::theme_bw()+theme(axis.title.y=element_blank()
                                                                       ,axis.text.y=element_text(family="serif",face="plain",size = 15)
                                                                       ,axis.text.x=element_text(family="serif",face="plain")
                                                                       ,plot.title = element_text(hjust = 0.5,family="serif",face="plain")
                                                                       ,axis.title.x=element_text(family="serif",face="plain")
                                                                       ,legend.title = element_text(family="serif",face="plain")
                                                                       ,legend.text = element_text(family="serif",face="plain"))
    }else if(colors=='nejm'){
      bubble=bubble2+ggsci::scale_fill_nejm()+ggplot2::theme_bw()+theme(axis.title.y=element_blank()
                                                                        ,axis.text.y=element_text(family="serif",face="plain",size = 15)
                                                                        ,axis.text.x=element_text(family="serif",face="plain")
                                                                        ,plot.title = element_text(hjust = 0.5,family="serif",face="plain")
                                                                        ,axis.title.x=element_text(family="serif",face="plain")
                                                                        ,legend.title = element_text(family="serif",face="plain")
                                                                        ,legend.text = element_text(family="serif",face="plain"))
    }else if(colors=='lancet'){
      bubble=bubble3+ggsci::scale_fill_lancet()+ggplot2::theme_bw()+theme(axis.title.y=element_blank()
                                                                          ,axis.text.y=element_text(family="serif",face="plain",size = 15)
                                                                          ,axis.text.x=element_text(family="serif",face="plain")
                                                                          ,plot.title = element_text(hjust = 0.5,family="serif",face="plain")
                                                                          ,axis.title.x=element_text(family="serif",face="plain")
                                                                          ,legend.title = element_text(family="serif",face="plain")
                                                                          ,legend.text = element_text(family="serif",face="plain"))
    }else if(colors=='jama'){
      bubble=bubble4+ggsci::scale_fill_jama()+ggplot2::theme_bw()+theme(axis.title.y=element_blank()
                                                                        ,axis.text.y=element_text(family="serif",face="plain",size = 15)
                                                                        ,axis.text.x=element_text(family="serif",face="plain")
                                                                        ,plot.title = element_text(hjust = 0.5,family="serif",face="plain")
                                                                        ,axis.title.x=element_text(family="serif",face="plain")
                                                                        ,legend.title = element_text(family="serif",face="plain")
                                                                        ,legend.text = element_text(family="serif",face="plain"))
    }else if(colors=='jco'){
      bubble=bubble5+ggsci::scale_fill_jco()+ggplot2::theme_bw()+theme(axis.title.y=element_blank()
                                                                       ,axis.text.y=element_text(family="serif",face="plain",size = 15)
                                                                       ,axis.text.x=element_text(family="serif",face="plain")
                                                                       ,plot.title = element_text(hjust = 0.5,family="serif",face="plain")
                                                                       ,axis.title.x=element_text(family="serif",face="plain")
                                                                       ,legend.title = element_text(family="serif",face="plain")
                                                                       ,legend.text = element_text(family="serif",face="plain"))
    }
    bubble=bubble+ggtitle(tl)
    bubble=bubble+theme(text=element_text(size=12,family="serif"))
  }
  return(bubble)
}

up.result=fun_clusterProfiler(up.genes,pvalueCutoff = 0.5,fdrCutoff = 0.5,qvalueCutoff = 0.5)
down.result=fun_clusterProfiler(down.genes,pvalueCutoff = 0.5,fdrCutoff = 0.5,qvalueCutoff = 0.5)

a=dotplot_batch(up.result$KEGG@result,db = 'KEGG pathway (Up)',colors = 'nejm')
b=dotplot_batch(up.result$GO_BP@result,db = 'GO (Up)',colors = 'nejm')
c=dotplot_batch(down.result$KEGG@result,db = 'KEGG pathway (Down)',colors = 'nejm')
d=dotplot_batch(down.result$GO_BP@result,db = 'GO (Down)',colors = 'nejm')
fuji=ggpubr::ggarrange(a,b,c,d, ncol = 2, nrow = 2,labels = '',align = "hv")

