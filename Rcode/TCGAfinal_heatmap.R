

library(data.table)
library(ComplexHeatmap)
library(dplyr)
options(stringsAsFactors = F)
DEGs<-fread("../dataset/TCGA_data/microSignature_DEGtag.csv")%>%as.data.frame()
exprall<-fread("../dataset/TCGA_data/expr_MicroSignature_TMP.csv")%>%as.data.frame()
sampleInfo<-fread("../dataset/TCGA_data/sampleInfo.csv")%>%as.data.frame()
MRscoreall<-fread("../dataset/TCGA_data/MRscore_TCGA_Patients9359.csv")%>%as.data.frame()
colnames(sampleInfo)[1]="sampleID"
colnames(MRscoreall)[1]="sampleID"
sampleInfo$sampleID<-gsub("-",".",sampleInfo$sampleID)
MRscoreall$sampleID<-gsub("-",".",MRscoreall$sampleID)
colnames(exprall)<-gsub("-",".",colnames(exprall))
DEGs=DEGs[,c(1,6,7)]
exprall[1:3,1:3]
sampleInfo[1:3,1:3]
MRscoreall[1:3,1:3]
head(DEGs)
cancerTypes="STAD"
load(file = "../dataset/MRscoreAPPexample.RData")
#exprMat=exprMat.exmaple
MRgenes=MRgenes.example
MRgeneAnno=MRgeneAnno.example

fun_to_heatmap<-function(cancerTypes){
  MRDEs<-subset(DEGs,Types==cancerTypes)
  meta<-subset(sampleInfo,Types==cancerTypes)
  expr<-exprall[which(exprall$Gene%in%MRDEs$Gene),]%>%as.data.frame()
  expr<-data.frame(row.names  = expr$Gene,expr[,which(colnames(expr)%in%meta$sampleID)])
  MRscore<-subset(MRscoreall,Types==cancerTypes)[,-2]
  expr[1:3,1:3]
  exprMat=expr
  MRgenes=MRgenes.example
  MRgeneAnno=MRgeneAnno.example
  antiVir<-MRgeneAnno[grep("vir",MRgeneAnno$GOterm),]
  antiBac<-MRgeneAnno[c(grep("bacter",MRgeneAnno$GOterm)
                        ,grep("lipopolysaccharide",MRgeneAnno$GOterm)),]
  MRgeneExpr<-expr
  antiVirResponse<-MRgeneExpr[which(rownames(MRgeneExpr)%in%antiVir$Gene),]
  antiBacResponse<-MRgeneExpr[which(rownames(MRgeneExpr)%in%antiBac$Gene),]
  MRscore<-data.frame(row.names = MRscore$sampleID,MRscore=MRscore$MRscore)
  dfcol<-data.frame(row.names =meta$sampleID,Condition=meta$Condition)
  annInfo=bind_cols(MRscore,antiVirResponse)
  ha = HeatmapAnnotation(
    Condition=anno_histogram(dfcol,Condition =c("Cancer" =  "red", 
                                                  "Normal" = "black")),
    antiVirResponse=anno_density(as.matrix(antiVirResponse), type =  "heatmap",
                                 gp = gpar(col = "white")),
    height = unit(2.5, "cm"),
    antiBacResponse=anno_density(as.matrix(antiBacResponse),type =  "heatmap")
  )
  
  p <-Heatmap(scale(log(MRgeneExpr+1),center = F),show_column_names = F,
              heatmap_legend_param = list(title = "Z-score",
                                          legend_direction = "vertical",
                                          labels_gp = gpar(fontsize = 9)),
              row_dend_width = unit(10, "mm"),
              column_dend_height = unit(5, "mm"),
              row_names_gp = gpar(fontsize = 2),
              column_names_gp = gpar(fontsize = 6),
              bottom_annotation = ha)
  
  p
  return(p)
}