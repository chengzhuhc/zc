library(maftools)
load('Data/TCGA_SKCM_maf.RData')
clinical1=clinical[unique(substr(colnames(data),1,12)),]
comid1=get_cli_data(as.character(clinType[1]),data,can,clinical1)
data.gene=as.numeric(data[,comid])
Barcode=unique(data.maf@data$Tumor_Sample_Barcode)
data.cli=data.frame('Tumor_Sample_Barcode'=Barcode[match(comid,substr(Barcode,1,15))]
                ,'Groups'=c(ifelse(data.gene>=median(data.gene),'High exp','Low exp')))
data.cli=data.cli[which(data.cli$Tumor_Sample_Barcode!='NA'),]
write.csv(data.cli,file = paste(result_path,'/',fileNamedata,'_mut_group.csv',sep = ""),quote = F,row.names = F)
data.maf@clinical.data$Tumor_Sample_Barcode=NULL
data.maf@clinical.data$Tumor_Sample_Barcode=data.cli$Tumor_Sample_Barcode
data.maf@clinical.data$Groups=data.cli$Groups
final.maf=data.maf
col = RColorBrewer::brewer.pal(n = 10, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Ins','In_Frame_Ins', 'Splice_Site', 'In_Frame_Del','Nonstop_Mutation','Translation_Start_Site','Multi_Hit')

plotmafSummary(maf=data.maf, rmOutlier=TRUE, addStat="median", dashboard=TRUE, titvRaw = FALSE,color = col)
oncoplot(maf = final.maf,genes = final.maf@gene.summary$Hugo_Symbol[1:20],clinicalFeatures = 'Groups',colors = col,sortByAnnotation = TRUE)
lollipopPlot(maf = final.maf, gene = dnaName, colors = col)


