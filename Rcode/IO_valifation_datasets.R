library(tidyverse)
library(data.table)
library(GEOquery)
options(stringsAsFactors = F)
GEOid<-c("GSE89997","GSE93157","GSE78220","GSE100797","GSE91061")
path1<-list.files("../dataset/IO_dataset/data/GSE/GSE89997_family/")
ges1<-path1[grep("tbl",path1)]
path2<-list.files("../dataset/IO_dataset/data/GSE/GSE93157_family/")
ges2<-path2[grep("tbl",path2)]
names1<-gsub("-tbl-1.txt","",ges1)
names2<-gsub("-tbl-1.txt","",ges2)
expr1<-list()
expr2<-list()
for (i in seq(length(ges1))) {
  expr1[[i]]<-fread(paste0("../dataset/IO_dataset/data/GSE/GSE89997_family/",ges1[i]))
  expr2[[i]]<-fread(paste0("../dataset/IO_dataset/data/GSE/GSE93157_family/",ges2[i]))
  names(expr1)[i]<-names1[i]
  names(expr2)[i]<-names2[i]
}
names(expr1)
names(expr2)
gpl1<-expr1[[1]]
expr1<-expr1[-1]
gpl2<-expr2[[1]]
expr2<-expr2[-1]
for (i in seq(length(expr1))) {
  colnames(expr1[[i]])<-c("probeID",names(expr1)[i])
}
for (i in seq(length(expr2))) {
  colnames(expr2[[i]])<-c("Gene",names(expr2)[i])
}
head(expr1[[1]])
head(expr2[[1]])
sapply(expr1,nrow)#70523 probes
sapply(expr2,nrow)#765 genes
expr1<-bind_cols(expr1)%>%as.data.frame()
expr1[1:5,1:5]
expr1<-expr1[,-grep("probeID[0-9]",colnames(expr1))]
test<-gpl1[1:5,1:15]
write.csv(test,"test.csv",row.names = F)
probeID<-gpl1[,c(1,8)]
colnames(probeID)<-c("probeID","ENSEMBLTRANS")
probeID$ENSEMBLTRANS<-gsub("^.*/// ENST","ENST",probeID$ENSEMBLTRANS)%>%gsub(" //.*","",.)
head(probeID)
new<-probeID[1:70753,1:2]
new<-new[grep("ENST",new$ENSEMBLTRANS),]
head(new)
library("clusterProfiler")
library(org.Hs.eg.db)
geneID <- bitr(new$ENSEMBLTRANS, fromType = "ENSEMBLTRANS", 
               toType =  c("ENSEMBLTRANS","SYMBOL"), 
               OrgDb = org.Hs.eg.db)
probe_anno<-merge(new,geneID,by="ENSEMBLTRANS")
head(probe_anno)
expr1[1:5,1:5]
expr1<-merge(probe_anno,expr1,by="probeID")
expr1<-expr1[,-c(1:2)]
colnames(expr1)[1]<-"Gene"
expr3<-merge(geneID,expr3,by="ENTREZID")
expr2<-bind_cols(expr2)%>%as.data.frame()
expr2[1:5,1:5]
expr2<-expr2[,-grep("Gene[0-9]",colnames(expr2))]
expr3<-fread("../dataset/IO_dataset/data/GSE/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv")
expr3[1:5,1:5]#gene is ENTREZID id need to transformed
colnames(expr3)[1]<-"ENTREZID"
library("clusterProfiler")
library(org.Hs.eg.db)
geneID <- bitr(expr3$Gene, fromType = "ENTREZID", 
                toType =  c("ENTREZID","SYMBOL"), 
                OrgDb = org.Hs.eg.db)
head(geneID)
expr3<-merge(geneID,expr3,by="ENTREZID")
expr3<-expr3[,-1]
colnames(expr3)[1]<-"Gene"
expr3[1:3,1:3]
expr4<-fread("../dataset/IO_dataset/data/GSE/GSE100797_ProcessedData.txt")
expr4[1:5,1:5]#checked
colnames(expr4)[1]<-"Gene"
expr5<-read.csv("../dataset/IO_dataset/data/GSE/GSE78220_PatientFPKM.csv")
expr5[1:5,1:5]#checked
expr4[1:5,1:5]
expr3[1:5,1:5]
expr2[1:5,1:5]
expr1[1:5,1:5]
#"GSE89997""GSE93157""GSE91061" "GSE100797" "GSE78220"
IO_dataset_exprMat<-list(GSE89997=expr1,#ARRAY
                         GSE93157=expr2,#ARRAY  
                         GSE91061=expr3,#FPKM
                         GSE100797=expr4,#ARRAY
                         GSE78220=expr5)#FPKM 
sapply(IO_dataset_exprMat, dim)
#GSE89997 GSE93157 GSE91061 GSE100797 GSE78220
#[1,]     5001      765    22086     18418    25268
#[2,]       31       31      110        26       29

meta2<-fread("../dataset/IO_dataset/data/GSE/GSE93157_family/SampleInfo_contrib1-GPL19965.txt")
meta5<-fread("../dataset/IO_dataset/data/GSE78220_series_matrix.txt")
meta4<-fread("../dataset/IO_dataset/data/GSE100797_series_matrix.txt")
      

