library(data.table)
library(tidyverse)
library(ggsci)
library(ggpubr)
library(ggthemes)
library(ggrepel)
options(stringsAsFactors = F)
files=list.files("../dataset/COX_results/")[grep("mutil_res_table_",list.files("../dataset/COX_results/"))]
files
files_single=list.files("../dataset/COX_results/")[grep("KM_res_single_TCGA_",list.files("../dataset/COX_results/"))]
files_single
types_single=gsub("KM_res_single_TCGA_","",files_single)%>%gsub(".csv","",.)
types_single
singleall<-list()
for (i in seq_along(files_single)) {
  singleall[[i]]<-fread(paste0(path,files_single[i]))%>%as.data.frame()%>%mutate(Types=rep(types_single[i]),nrow(.))
  names(singleall)[i]=types_single[i]
  colnames(singleall[[i]])[1]="Gene"
  
}
head(singleall[[1]])
km_MRgene<-bind_rows(singleall)
colnames(km_MRgene)[6]="geneNum"
head(km_MRgene)
coxall<-list()
path="../dataset/COX_results/"
for (i in seq_along(files)) {
  coxall[[i]]<-fread(paste0(path,files[i]))%>%as.data.frame()%>%mutate(Types=rep(types[i]),nrow(.))
  names(coxall)[i]=types[i]
  colnames(coxall[[i]])[1]="Gene"
   }
coxall[[1]]
sapply(coxall, nrow)
#ACC BLCA BRCA CESC CHOL COAD DLBC* ESCA  GBM* KICH KIRC* KIRP LAML LUAD LUSC* MESO   OV PRAD 
#8    8  122    1    4   16   29  106  133   13   13    2    5   61    3    3    6    1
sapply(singleall, nrow)
#ACC BLCA BRCA CESC CHOL COAD DLBC ESCA  GBM HNSC KICH KIRC KIRP LAML  LGG LIHC LUAD LUSC 
#119   92  128   12   67  161  283  117  149  298  203  128  138   50    6   10  214   55 
#MESO   OV PRAD SARC SKCM STAD 
#27  177    1  121   21    6 
MRgeneFC<-fread("../dataset/TCGA_data/microSignature_DEGtag.csv")
head(MRgeneFC)
length(unique(factor(km_MRgene$Gene)))
length(unique(factor(MRgeneFC$Gene)))
MRgeneKM<-merge(MRgeneFC,km_MRgene,by=c("Gene","Types"))
head(MRgeneKM)
MRgeneKM$HR<-MRgeneKM$`HR (95% CI for HR)`
MRgeneKM$HR<-gsub(" .*","",MRgeneKM$HR)%>%as.numeric()
head(MRgeneKM)
summary(MRgeneKM)
MRgeneKM$group<-MRgeneKM$HR
MRgeneKM$group<-ifelse(MRgeneKM$Log2FC<=-2&MRgeneKM$HR>=1.5,"HR1.5_FCd2",
                       ifelse(MRgeneKM$Log2FC<=-2&MRgeneKM$HR<=0.5,"HR0.5_FCd2",
                              ifelse(MRgeneKM$Log2FC>=2&MRgeneKM$HR>=1.5,"HR1.5_FCu2",   
                                     ifelse(MRgeneKM$Log2FC>=2&MRgeneKM$HR<=0.5,"HR0.5_FCu2","NO"))))
levels(factor(MRgeneKM$group))

p<-ggplot(MRgeneKM,aes(Log2FC,HR,color=group))+
  geom_point(aes(alpha = 0.6,size=-log10(p.value)), show.legend = FALSE)+
  scale_size(range = c(0, 4))+
  geom_vline(xintercept = c(-2, 2), color = 'gray', linetype = 2, size = 1) + 
  geom_hline(yintercept = c(0.5,1.5), color = 'gray', linetype = 2, size = 1)+
  scale_color_manual(limits = c('HR1.5_FCd2', 'HR0.5_FCd2', 'HR1.5_FCu2', 'HR0.5_FCu2', 'NO'), 
                     values = c('red', 'orange', 'purple', 'blue', 'gray'))+
  theme_classic2()
MRgeneKM$label=ifelse(abs(MRgeneKM$Log2FC)>=2&(MRgeneKM$HR>=1.5|MRgeneKM$HR<=0.5),MRgeneKM$Gene,"")

p+geom_text_repel(data=MRgeneKM,aes(Log2FC,HR,label=label),
                  size = 3,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE)
