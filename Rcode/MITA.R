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
genes<-levels(factor(keypathwaygeneExpr$Gene))##83genes
data1<-data1[which(data1$Gene%in%genes),]
mat<-data.frame(row.names = data1$Gene,data1[,-c(1:2)])
dfrow<-data.frame(row.names = data1$Gene,KEGG_term=data1$KEGG_term)
pheatmap(mat,border_color = NA,
         color = palette(value = c("blue","white","red")(100)),
         annotation_row = dfrow,
         cluster_rows = T)
