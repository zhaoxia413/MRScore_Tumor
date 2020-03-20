library(data.table)
library(ggsci)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
options(stringsAsFactors = F)
df<-fread("../dataset/TCGA_data/MRgenes_CancerTYpes_DEGtag.csv")
up<-subset(df,Regulation=="Up")
down<-subset(df,Regulation=="Down")
genelist_up<-split(up,f=up$Types)
genelist_down<-split(down,f=down$Types)
keypath<-fread("../dataset/TCGA_data/Nod_TLR_genes_1.csv")[,c(2,4)]
keypathwaygeneExpr<-merge(keypath,df,by="Gene")
upkey<-list()
downkey<-list()
for (i in seq_along(genelist_down)) {
  upkey[[i]]<-merge(keypath,genelist_up[[i]],by="Gene")
  downkey[[i]]<-merge(keypath,genelist_down[[i]],by="Gene")
  names(upkey)[i]<-names(genelist_up)[i]
  names(downkey)[i]<-names(genelist_down)[i]
}
head(upkey[[3]])
data<-keypathwaygeneExpr[,c(1,3,5)]
data$Log2FC<-round(data$Log2FC,2)
head(data)
datalist<-split.data.frame(data,f=data$Types,drop = F)
for (i in seq_along(datalist)) {
  colnames(datalist[[i]])[2]<-names(datalist)[i]
  datalist[[i]]<-datalist[[i]][,-3]
}
head(datalist[[1]])
head(keypath)
data1<-left_join(keypath,datalist[[1]],by="Gene",all=T)%>%
  left_join(.,datalist[[2]],by="Gene",all=T)%>%
  left_join(.,datalist[[3]],by="Gene",all=T)%>%
  left_join(.,datalist[[4]],by="Gene",all=T)%>%
  left_join(.,datalist[[5]],by="Gene",all=T)%>%
  left_join(.,datalist[[6]],by="Gene",all=T)%>%
  left_join(.,datalist[[7]],by="Gene",all=T)%>%
  left_join(.,datalist[[8]],by="Gene",all=T)%>%
  left_join(.,datalist[[9]],by="Gene",all=T)%>%
  left_join(.,datalist[[10]],by="Gene",all=T)%>%
  left_join(.,datalist[[11]],by="Gene",all=T)%>%
  left_join(.,datalist[[12]],by="Gene",all=T)%>%
  left_join(.,datalist[[13]],by="Gene",all=T)%>%
  left_join(.,datalist[[14]],by="Gene",all=T)%>%
  left_join(.,datalist[[15]],by="Gene",all=T)%>%
  left_join(.,datalist[[16]],by="Gene",all=T)%>%
  left_join(.,datalist[[17]],by="Gene",all=T)%>%
  left_join(.,datalist[[18]],by="Gene",all=T)%>%
  left_join(.,datalist[[19]],by="Gene",all=T)%>%
  left_join(.,datalist[[20]],by="Gene",all=T)%>%
  left_join(.,datalist[[21]],by="Gene",all=T)%>%
  left_join(.,datalist[[22]],by="Gene",all=T)%>%
  left_join(.,datalist[[23]],by="Gene",all=T)%>%
  left_join(.,datalist[[24]],by="Gene",all=T)%>%
  left_join(.,datalist[[25]],by="Gene",all=T)%>%
  left_join(.,datalist[[26]],by="Gene",all=T)%>%
  left_join(.,datalist[[27]],by="Gene",all=T)%>%
  left_join(.,datalist[[28]],by="Gene",all=T)
data1[is.na(data1)]=0
data1<-data1[c(1:231),]
test<-data1%>%as.data.frame()%>%mutate(sum=sum())
library(pheatmap)
library(RColorBrewer)
genes<-levels(factor(keypathwaygeneExpr$Gene))##83genes
data1<-data1[which(data1$Gene%in%genes),]
write.csv(data1,"MRgenes_TCGA_keygenes.csv",row.names = F)
mat<-data.frame(row.names = data1$Gene,data1[,-c(1:2)])
dfcol<-data.frame(row.names = data1$Gene,KEGG_term=data1$KEGG_term)
bk <- c(seq(-9.39,0,by=0.1),seq(0,6.84,by=0.1))
pheatmap(mat,border_color = NA,
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","orange","red"))(length(bk)/2)),
         annotation_row = dfrow,
         cluster_rows = T)
library(ConsensusClusterPlus)
dir.create("ConsensusClusterPlusresults")
title="./ConsensusClusterPlusresults"
d<-as.matrix(mat)
#d<-log(d+1)
d = sweep(d,1, apply(d,1,median,na.rm=T))
results<-ConsensusClusterPlus(d,maxK=8,reps=50,
                              pItem=0.8,pFeature=1,
                              title=title,clusterAlg="hc",
                              distance="spearman",
                              seed=1262118388.71279,plot="png")
cluster<-results[[4]]["consensusClass"]
write.csv(cluster,"ConsensusCluster_byCancers.csv")

color.key=c("#112D4E","#0245A3","#8FBAF3","#BDF1F6","white","#FF6F3C","#E84545","#D72323","#D72323")


clus<-read.csv("ConsensusCluster_byCancers.csv",header = T)%>%as.data.frame()
data2<-t(mat)
data2<-data.frame(Types=rownames(data2),data2)
data3<-merge(clus,data2,by="Types")
data3<-data3[order(data3$consensusClass),]
mat<-data.frame(row.names = data3$Types,data3[,-c(1:2)])
dfrow<-data.frame(row.names = data3$Types,consensusClass=data3$consensusClass)
ann_colors = list(
  consensusClass = c(Cluster1 = "gray", Cluster2 ="orange",
                     Cluster3 ="darkgreen",Cluster4 = "#252A34"),
  KEGG_term=c("NLR/TLR" = "#FFF2BE",
              NLR= "#F6318C",
              TLR="blue"))
p<-pheatmap(mat,
            #kmeans_k = 8,
            color = colorRampPalette(color.key)(15),
            border_color = NA,
            annotation_row = dfrow,
            annotation_col = dfcol,
            cluster_rows = F,
            cluster_cols = T,
            #clustering_distance_cols = "euclidean",
            scale="row",
            #clustering_distance_rows="euclidean",
            annotation_colors = ann_colors,
            #clustering_method = "complete",
            show_rownames = T,
            show_colnames = T,
            annotation_names_row = T,
            annotation_names_col = T)
##micro_immunedatasets
files<-list.files("../dataset/dataset_alidation/validation_results/")
event<-files[grep("csv",files)]
event<-event[-c(2:4,6)]
event
filenames<-gsub("_DEGres.csv","",event)
DEGs<-list()
sigDEGs<-list()
for (i in seq_along(event)) {
  DEGs[[i]]<-fread(paste0("../dataset/dataset_alidation/validation_results/",event[i]))%>%as.data.frame()
  sigDEGs[[i]]<-subset(DEGs[[i]],adj.P.Val<=0.05)
  names(sigDEGs)[i]<-filenames[i]
  colnames(sigDEGs[[i]])[c(1,2)]<-c("Gene",filenames[i])
  sigDEGs[[i]]<-sigDEGs[[i]][c(1,2)]
}
sigDEGs[[1]][1:3,]
keypath<-fread("../dataset/TCGA_data/Nod_TLR_genes_1.csv")[,c(2,4)]
data1<-left_join(keypath,sigDEGs[[1]],by="Gene",all=T)%>%
  left_join(.,sigDEGs[[2]],by="Gene",all=T)%>%
  left_join(.,sigDEGs[[3]],by="Gene",all=T)%>%
  left_join(.,sigDEGs[[4]],by="Gene",all=T)%>%
  left_join(.,sigDEGs[[5]],by="Gene",all=T)%>%
  left_join(.,sigDEGs[[6]],by="Gene",all=T)
data1[is.na(data1)]=0
data2<-data1%>%mutate(sum=apply(data1[,-c(1,2)], 1, sum))
data2<-subset(data2,sum>0)##63 genes
write.csv(data2,"MRgenes_microbialinfection_keygenes.csv",row.names = F)
mat<-data.frame(row.names = data2$Gene,round(data2[,-c(1:2,9)],2))
p<-pheatmap(mat,
            scale = "none",
            color = colorRampPalette(color.key)(15),
            border_color = NA)
##LUAD datasets
load(file = "../dataset/dataset_alidation/LUSC/LUSC_9datasets_expr_meta.Rdata")
expr<-LUSC_9datasets_expr_meta$LUSC_expr
meta<-LUSC_9datasets_expr_meta$LUSC_meta

meta[1:3,1:10]
DEG<-fread("../dataset/LUAD_further/DEGs.csv")%>%as.data.frame()
head(DEG)
colnames(DEG)[1]<-"Gene"
DEGs_filter<-subset(DEG,abs(logFC)>=1&P.Value<=0.05)
expr[1:3,1:3]
DEGs_filter[1:3,1:3]
DEexpr<-DEGs_filter
expr2MRscore<-function(expr,DEexpr){
  signature<-fread("../dataset/TCGA_data/annotationRow1.csv")
  MRgene_DEGs<-subset(DEexpr,Gene%in%signature$Gene)
  MRgene_DEGs$Regulation<-sample(c("Up","Down"),nrow(MRgene_DEGs),replace = T)
  MRgene_DEGs$Regulation<-ifelse(MRgene_DEGs$logFC>=0,"Up","Down")
  DE_Mgene<-MRgene_DEGs[,c("Gene","Regulation")]
  Total_MRscore<-sum(subset(MRgene_DEGs,logFC>=0)$logFC)+sum(subset(MRgene_DEGs,logFC<=0)$logFC)
  message(paste0("Total Score = ",Total_MRscore))
  upgene<-subset(DE_Mgene,Regulation=="Up")
  downgene<-subset(DE_Mgene,Regulation=="Down")
  Upexpr<-expr[which(expr$Gene%in%upgene$Gene),]
  Downexpr<-expr[which(expr$Gene%in%downgene$Gene),]
  test<-Upexpr[,-grep("Gene",colnames(Upexpr))]
  test
  apply(test,2,mean)
  Upexpr_average<-data.frame(sampleID = colnames(Upexpr)[-grep("Gene",colnames(Upexpr))],Upscore=apply(Upexpr[,-grep("Gene",colnames(Upexpr))], 2, mean))
  Downexpr_average<-data.frame(sampleID = colnames(Downexpr)[grep("Gene",colnames(Downexpr))],Downscore=apply(Downexpr[,-grep("Gene",colnames(Downexpr))], 2, mean))
  DE_average<-bind_cols(Upexpr_average,Downexpr_average)[,-3]
  MRscore_value<-DE_average%>%mutate(MRscore=Upscore-Downscore)
  return(list(MRgene_DEGs=MRgene_DEGs,MRscore_value=MRscore_value,Total_MRscore=Total_MRscore))
}  
MRscorelist<-expr2MRscore(expr = expr,DEexpr = DEGs_filter)  
MRscore<-MRscorelist[[2]]
write.csv(MRscorelist[[1]],"MIRenes_LUAD.csv",row.names = F)
head(MRscore)
