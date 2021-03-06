library(data.table)
library(tidyverse)
library(ggsci)
library(ggpubr)
library(ggthemes)
library(ggrepel)
options(stringsAsFactors = F)
files=list.files("../dataset/COX_results/")[grep("mutil_res_table_",list.files("../dataset/COX_results/"))]
files
path="../dataset/COX_results/"
types=gsub("mutil_res_table_","",files)%>%gsub(".csv","",.)
types
files_single=list.files("../dataset/COX_results/")[grep("KM_res_single_TCGA_",list.files("../dataset/COX_results/"))]
files_single
types_single=gsub("KM_res_single_TCGA_","",files_single)%>%gsub(".csv","",.)
types_single

coxall<-list()
for (i in seq_along(files)) {
  print(i)
  coxall[[i]]<-fread(paste0(path,files[i]))%>%as.data.frame()%>%mutate(Types=rep(types[i]),nrow(.))
  names(coxall)[i]=types[i]
  colnames(coxall[[i]])[c(1,3,4)]=c("Gene","HR","se")
}
singleall<-list()
for (i in seq_along(files)) {
  print(i)
  singleall[[i]]<-fread(paste0(path,files_single[i]))%>%as.data.frame()%>%mutate(Types=rep(types_single[i]),nrow(.))
  names(singleall)[i]=types_single[i]
  colnames(singleall[[i]])[1]="Gene"
}

sapply(coxall, nrow)
coxall<-coxall[-which(sapply(coxall, nrow)==1)]
cox_MRgene<-bind_rows(coxall[c(1:6,9:16)])
#ACC BLCA BRCA CESC CHOL COAD DLBC* ESCA  GBM* KICH KIRC* KIRP LAML LUAD LUSC* MESO   OV PRAD 
#8    8  122    1    4   16   29  106  133   13   13    2    5   61    3    3    6    1
sapply(singleall, nrow)
#ACC BLCA BRCA CESC CHOL COAD DLBC ESCA  GBM HNSC KICH KIRC KIRP LAML  LGG LIHC LUAD LUSC 
#119   92  128   12   67  161  283  117  149  298  203  128  138   50    6   10  214   55 
#MESO   OV PRAD SARC SKCM STAD 
#27  177    1  121   21    6 
MRgeneFC<-fread("../dataset/TCGA_data/microSignature_DEGtag.csv")
head(MRgeneFC)
km_MRgene<-bind_rows(singleall)
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
keygenes=subset(MRgeneKM,label!="")
colnames(keygenes)
<<<<<<< HEAD
p=ggdotchart(keygenes, x="Gene", y="HR", color = "Types",
           #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
=======
keygenes<-keygenes[order(keygenes$Gene,keygenes$HR)]

p=ggdotchart(keygenes, x="Gene", y="coef", color = "Types",fill="Types",
           palette = col11,
>>>>>>> 0f9b27534ad7178c5fb10cb34e46741cc72c6ced
           sorting = "descending", add = "segments", rotate = TRUE,
           group = "Types", dot.size = 6,shape = factor(keygenes$Log2FC>=0),
           label = round(keygenes$HR,2), font.label = list(color="black",
                                                     size=9, vjust=0.5), ggtheme = theme_pubr())
p+geom_hline(yintercept = 0, linetype = 2, color = "lightgray")
library(clusterProfiler)
library(org.Hs.eg.db)
upgenes<-subset(keygenes,Regulation>=0)$Gene%>%as.factor()
downgenes<-subset(keygenes,Regulation<=0)$Gene%>%as.factor()
genes_ID=select(org.Hs.eg.db,keys =levels(factor(keygenes$Gene)),column = "ENTREZID", keytype = "SYMBOL",multiVals = "first")
genes_ID=genes_ID[,2]
ekk <- enrichKEGG(gene = genes_ID,organism = 'hsa',pvalueCutoff = 0.05)
write.table(as.matrix(ekk@result), 
            file=paste0("../dataset/COX_results/KEGGenrich.txt"),
            row.names = F,sep = "\t",quote = F)
df<-ekk@result%>%filter(qvalue<=0.01)%>%mutate(Bg=gsub("/.*","",BgRatio))%>%
  mutate(enrichScore = (Count/as.numeric(Bg)))%>%  mutate(NegLog_qvalue = -log10(qvalue))
colnames(df)
df<-df[order(df$enrichScore),]
p<-ggplot(df,aes(enrichScore,reorder(Description,enrichScore)))+
  geom_point(aes(size=Count,color=NegLog_qvalue))+ 
  scale_color_gradient(low = "blue",high = "red")+
  theme(axis.text = element_text(colour = "black",size = 12),
        legend.text = element_text(colour = "black",size = 12),
        panel.grid.major = element_line(colour = "gray88",linetype = "dashed"),
        panel.background=element_blank(),
        axis.ticks = element_line(size = 1, colour = "black"),
        axis.line = element_line(size = 1, colour = "black"),
        legend.position = "right",
        panel.border = element_rect(linetype = "solid", fill = NA,size = 2),
        axis.title = element_text(size = 12)
        ,title = element_text(size = 12))
p
