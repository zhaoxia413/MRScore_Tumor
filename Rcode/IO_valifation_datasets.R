library(tidyverse)
library(data.table)
library(GEOquery)
options(stringsAsFactors = F)
GEOid<-c("GSE93157","GSE78220","GSE100797","GSE91061")
path2<-list.files("../dataset/IO_dataset/data/GSE/GSE93157_family/")
ges2<-path2[grep("tbl",path2)]
names2<-gsub("-tbl-1.txt","",ges2)
expr2<-list()
for (i in seq(length(ges2))) {
  expr2[[i]]<-fread(paste0("../dataset/IO_dataset/data/GSE/GSE93157_family/",ges2[i]))
  names(expr2)[i]<-names2[i]
}
names(expr2)
gpl2<-expr2[[1]]
expr2<-expr2[-1]
for (i in seq(length(expr2))) {
  colnames(expr2[[i]])<-c("Gene",names(expr2)[i])
}
head(expr2[[1]])
sapply(expr2,nrow)#765 genes
expr2<-bind_cols(expr2)%>%as.data.frame()
expr2[1:5,1:5]
expr2<-expr2[,-grep("Gene[0-9]",colnames(expr2))]
expr3<-fread("../dataset/IO_dataset/data/GSE/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv")
expr3[1:5,1:5]#gene is ENTREZID id need to transformed
colnames(expr3)[1]<-"ENTREZID"
library("clusterProfiler")
library(org.Hs.eg.db)
geneID <- bitr(expr3$ENTREZID, fromType = "ENTREZID", 
                toType =  c("ENTREZID","SYMBOL"), 
                OrgDb = org.Hs.eg.db)
head(geneID)
expr3<-merge(geneID,expr3,by="ENTREZID")
expr3<-expr3[,-1]
colnames(expr3)[1]<-"Gene"
expr3[1:3,1:3]
expr4<-fread("../dataset/IO_dataset/data/GSE100797_ProcessedData.txt")
expr4[1:5,1:5]#checked
colnames(expr4)[1]<-"Gene"
expr5<-read.csv("../dataset/IO_dataset/data/GSE/GSE78220_PatientFPKM.csv")
expr5[1:5,1:5]#checked
expr4[1:5,1:5]
expr3[1:5,1:5]
expr2[1:5,1:5]
#"GSE93157""GSE91061" "GSE100797" "GSE78220"
IO_dataset_exprMat<-list(GSE93157=expr2,#ARRAY  
                         GSE91061=expr3,#FPKM
                         GSE100797=expr4,#ARRAY
                         GSE78220=expr5)#FPKM 
sapply(IO_dataset_exprMat, dim)
# GSE93157 GSE91061 GSE100797 GSE78220
#[1,] 765    22086     18418    25268
#[2,] 31       31      110        26       29
meta2<-fread("../dataset/IO_dataset/data/GSE/GSE93157_family/SampleInfo_contrib1-GPL19965.txt")
meta3<-fread("../dataset/IO_dataset/data/GSE91061_meta.txt")
meta5<-fread("../dataset/IO_dataset/data/GSE78220_meta.txt")
gset<-getGEO('GSE100797',destdir='.',
             AnnotGPL=F,
             getGPL=F)
a=gset[[1]] 
class(a)
dat=exprs(a) #a现在是一个对象，取a这个对象通过看说明书知道要用exprs这个函数
meta4<-pData(a) #用pData来提取临床信息
meta2[1:3,1:3]
meta3[1:3,1:3]
meta4[1:3,1:3]
meta5[1:3,1:3]
IO_dataset_meta<-list(GSE93157=meta2,#ARRAY  
           GSE91061=meta3,#FPKM
           GSE100797=meta4,#ARRAY
           GSE78220=meta5)#FPKM 
group2<-data.frame(sampleID=meta2$Accession,
                   cancerType=meta2$`1:Source`,
                   Treatment=meta2$`1:DRUG`,
                   PFS=meta2$`1:PFS`,
                   Status_PFS=meta2$`1:PFSE`,
                   Response=meta2$`1:RESPONSE`,
                   Best_Response=meta2$`1:BEST.RESP`)
group3<-data.frame(sampleID=meta3$acc,
                   Response=meta3$response,
                   cancerType=meta3$tissue)
group4<-data.frame(sampleID=meta4$geo_accession,
                   Response=meta4$characteristics_ch1.9,
                   TumorSize=meta4$`tumor.response:ch1`,
                   Treatment=meta4$characteristics_ch1,
                   PFS=meta4$`pfs.time:ch1`,
                   Status_PFS=meta4$`pfs.event:ch1`,
                   OS=meta4$`os.time:ch1`,
                   Status_OS<-meta4$`os.event:ch1`,
                   cancerType=meta4$source_name_ch1)
group5<-data.frame(sampleID=meta5$Run,
                   Treatment=meta5$treatment,
                   Response=meta5$`anti-pd-1_response`,
                   Status_OS=meta5$vital_status)
head(group5)
levels(factor(group5$Response))
group5$Response<-ifelse(group5$Response=="Complete Response","CR",
                        ifelse(group5$Response=="Partial Response","PR","PD"))
IO_dataset_group<-list(GSE93157=group2,#ARRAY  
                            GSE91061=group3,#FPKM
                            GSE100797=group4,#ARRAY
                            GSE78220=group5) #FPKM )
IO_dataset<-list(IO_dataset_exprMat,IO_dataset_meta,IO_dataset_group)
save(IO_dataset,file = "../dataset/IO_dataset/data/IO_dataset_expr_meta_group_list.Rdata" )
load(file = "../dataset/IO_dataset/data/IO_dataset_expr_meta_group_list.Rdata")
exprlist<-IO_dataset$exprMat
grouplist<-IO_dataset$group
names(exprlist)
#[1] "GSE93157"  "GSE91061"  "GSE100797" "GSE78220" 
#####samples
# ARRAY     FPKM     ARRAY      FPKM
# GSE93157 GSE91061 GSE100797 GSE78220
#[1,] 765    22086     18418    25268
head(grouplist[[1]])##pfs  ARRAY Melanoma
head(grouplist[[2]])##NO PFS or OS,response FPKM Melanoma
head(grouplist[[3]])#PFS and OS,response ARRAY Melanoma
head(grouplist[[4]])##only Status_OS,response FPKM Melanoma



