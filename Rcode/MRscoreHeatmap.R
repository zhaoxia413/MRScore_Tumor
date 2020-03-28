library(pheatmap)
library(data.table)
library(tidyverse)
library(ggsci)
library(ggthemes)
library(ggpubr)
library(RColorBrewer)
col1<-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                        "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
col31<-c("#303841","#D72323","#377F5B","#375B7F","#F2FCFC","#f0027f",
         "#FAF8DE","#666666","#BDF1F6","#023782","#5e4fa2","#F1C40F",
         "#ff7f00","#cab2d6","#240041","#ffff99","#0E3BF0","#a65628",
         "#f781bf","#808FA6","#2EB872","#F0FFE1","#F33535","#011F4E",
         "#82B269","#D3C13E","#3F9DCD","#014E1F","#AFFFDF","#3D002E",
         "#3A554A")
MRscore<-fread("../dataset/TCGA_data/MRscore_TCGA_Patients9359.csv")%>%as.data.frame()
ggboxplot(MRscore,x="Types",y="MRscore",fill="Types")+
  theme_few(base_size = 12)+
    geom_jitter(alpha=0.1,size=0.5)+
  theme(axis.title.y = element_blank())+
  scale_fill_manual(values = col31)+
  guides(fill=F)+
  coord_flip()
MRscore$Types<-MRscore$Types[order(MRscore$MRscore)]
MRscore$Types<-MRscore$Types[order(MRscore$Types,decreasing = T)]
ggboxplot(MRscore,x="Types",y="MRscore",fill="Types")+
  theme_few(base_size = 12)+
  geom_jitter(alpha=0.05,size=1)+
  theme(axis.title.y = element_blank())+
  scale_fill_manual(values = col31)+
  guides(fill=F)+
  coord_flip()
library(reshape2)
head(MRscore)
data<-dcast(MRscore,sample~Types)
data[is.na(data)]=0
mat<-data.frame(row.names = data$sample,as.matrix(data[,-1]))
pheatmap(scale(t(mat)),show_colnames  = F)
  


