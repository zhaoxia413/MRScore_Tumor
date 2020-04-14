

library(data.table)
library(ComplexHeatmap)
library(dplyr)
options(stringsAsFactors = F)
DEGs<-fread("../dataset/TCGA_data/microSignature_DEGtag.csv")%>%
  as.data.frame()
DEGs%>%group_by(Types)%>%summarise(n())

exprall<-fread("../dataset/TCGA_data/expr_MicroSignature_TMP.csv")%>%as.data.frame()
sampleInfo<-fread("../dataset/TCGA_data/sampleInfo.csv")%>%as.data.frame()
MRscoreall<-fread("../dataset/TCGA_data/MRscore_TCGA_Patients9359.csv")%>%as.data.frame()
colnames(sampleInfo)[1]="sampleID"
colnames(MRscoreall)[1]="sampleID"
sampleInfo$sampleID<-gsub("-",".",sampleInfo$sampleID)
MRscoreall$sampleID<-gsub("-",".",MRscoreall$sampleID)
colnames(exprall)<-gsub("-",".",colnames(exprall))
DEGs=DEGs[,c(1,4,6,7)]
exprall[1:3,1:3]
sampleInfo[1:3,1:3]
MRscoreall[1:3,1:3]
head(DEGs)
cancerTypes="BRCA"
library(data.table)
library(ComplexHeatmap)
library(dplyr)
load(file = "../dataset/MRscoreAPPexample.RData")
#clinic<-fread("../dataset/TCGA_data/TCGA_clinical.csv")%>%as.data.frame()
#colnames(clinic)[1]<-"sampleID"
#clinic$sampleID<-gsub("-",".",clinic$sampleID)
#clinic1<-fread("../dataset/TCGA_data/TCGA_Clinical_Variables.txt")%>%as.data.frame()
#clinic2<-clinic1[c(162,7,48,50,56,22,24),]
#rownames(clinic2)<-clinic2$V1
#clinic2<-data.frame(t(clinic2))[-1,]
#clinic2$sampleID<-rownames(clinic2)
#clinic2$sampleID<-gsub("-",".",clinic2$sampleID)
#clinic2$sampleID[1:3]
#clinic7706=clinic2
#clinic11160=clinic
#colnames(clinic7706)=c("Age","Gender","Tumor_pt","Nodes_pn","Metastasis_pm","OS_status","OStime","sampleID")
#save(DEGs,clinic7706,clinic11160,exprall,sampleInfo,MRscoreall,file = "../dataset/TCGA_results/TCGAMRegnesHeatmap.Rdata")
load(file = "../dataset/TCGA_results/TCGAMRegnesHeatmap.Rdata")
MRgenes<-MRgenes.example
head(sampleInfo,3)
MRscoreall[1:3,1:3]
MRgenes[1:3,1:4]
exprall[1:3,1:3]
clinic11160$sampleID[1:3]
clinic11160$sampleID1=clinic11160$sampleID
sampleInfo$sampleID1=str_sub(sampleInfo$sampleID,1,12)
meta11160<-merge(sampleInfo,clinic11160,by="sampleID1")
exprall1<-exprall
colnames(exprall1)=str_sub(colnames(exprall1),1,12)
clinic7706$sampleID[1:3]
meta7706<-merge(sampleInfo,clinic7706,by="sampleID")
cancer="STAD"
MRgenes_1<-subset(MRgenes,Types==cancer)
meta_1<-subset(meta7706,Types==cancer)
expr_1<-data.frame(row.names = exprall$Gene,exprall[,colnames(exprall)%in%meta_1$sampleID])
epr_1<-data.frame(t(expr_1))
epr_1<-data.frame(sampleID=rownames(epr_1),epr_1)
tcgaMLdt_1<-merge(meta_1,epr_1,by="sampleID")


meta1_1<-subset(meta11160,Types==cancer)
expr1_1<-data.frame(row.names = exprall1$Gene,exprall1[,colnames(exprall1)%in%meta1_1$sampleID1])
epr1_1<-data.frame(t(expr1_1))
epr1_1<-data.frame(sampleID1=rownames(epr1_1),epr1_1)
tcgaMLdt1_1<-merge(meta1_1,epr1_1,by="sampleID1")
levels(factor(tcgaMLdt1_1$Condition))
write.csv(tcgaMLdt1_1,"../dataset/MRgene_selection/tcgaMLdt_STAD.csv",row.names = F)



#exprMat=exprMat.exmaple
MRgenes=MRgenes.example
MRgeneAnno=MRgeneAnno.example
nsample=50
fun_to_heatmap("BRCA",FC = 4)
FC=2
cancerTypes="BRCA"
fun_to_heatmap<-function(cancerTypes,nsample,FC){
  DEGs<-filter(DEGs,nsample,abs(Log2FC)>=FC)
  MRDEs<-subset(DEGs,Types==cancerTypes)
  meta<-subset(sampleInfo,Types==cancerTypes)
  meta<-meta[sample(c(1:nrow(meta)),nsample),]
  expr<-exprall[which(exprall$Gene%in%MRDEs$Gene),]%>%as.data.frame()
  expr<-data.frame(row.names  = expr$Gene,expr[,which(colnames(expr)%in%meta$sampleID)])
  MRscore<-subset(MRscoreall,Types==cancerTypes)[,-2]
  exprMat=expr
  MRgenes=MRgenes.example
  MRgeneAnno=MRgeneAnno.example
  antiVir<-MRgeneAnno[grep("vir",MRgeneAnno$GOterm),]
  antiBac<-MRgeneAnno[c(grep("bacter",MRgeneAnno$GOterm)
                        ,grep("lipopolysaccharide",MRgeneAnno$GOterm)),]
  MRgeneExpr<-exprMat
  antiVirResponse<-MRgeneExpr[which(rownames(MRgeneExpr)%in%antiVir$Gene),]
  antiBacResponse<-MRgeneExpr[which(rownames(MRgeneExpr)%in%antiBac$Gene),]
  MRscore<-data.frame(row.names = MRscore$sampleID,MRscore=MRscore$MRscore)
  meta$sampleID=meta$sampleID[order(meta$Condition)]
  dfcol<-data.frame(meta$sampleID[order(meta$Condition)],
                    meta$Condition[order(meta$Condition)])
  dfcol1<-data.frame(dfcol[,2])
  rownames(dfcol1) <- dfcol$meta.sampleID.order.meta.Condition..
  colnames(dfcol1)="Condition"
  ha1= HeatmapAnnotation(df=dfcol1, which = "column",
                         col = list(Condition = c("Cancer" =  "orange", 
                                        "Normal" = "black")))
  ha = HeatmapAnnotation(
    antiVirResponse=anno_density(as.matrix(antiVirResponse), type =  "heatmap",
                                 gp = gpar(col = "white")),
    height = unit(2.5, "cm"),
    antiBacResponse=anno_density(as.matrix(antiBacResponse),type =  "heatmap")
  )
  
  p <-Heatmap(scale(log(MRgeneExpr+1),center = F),
              top_annotation = ha1,
              cluster_columns = F,
              show_column_names = F,
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
cancers<-levels(factor(DEGs$Types))
p<-list()
for (i in c(2:18,21:23,27:28)) {
  print(i)
  p[[i]]<-fun_to_heatmap(cancerTypes = cancers[i],nsample = 10,FC = 1)
names(p)[i]<-cancers[i]
  }

fun_to_heatmap(cancerTypes = "ACC",nsample = 10,FC = 1)
p[[3]]
p[[2]]
p[[4]]
p[[5]]
  