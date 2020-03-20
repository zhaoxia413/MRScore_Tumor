library(data.table)
library(ggsci)
library(ggplot2)
library(ggthemes)
  library(ggrepel)
OR<-fread("../dataset/OR/ORR29.csv")%>%as.data.frame()
MIR<-fread("../dataset/OR/Microscores_Cancertypes.csv")%>%as.data.frame()
data<-merge(OR,MIR,by="Types")
head(data)
data1<-data[-which(is.na(data$OS)),]
ggplot(data,aes(MRscore,ORR))+
  geom_point(aes(size=N,color=OS))+
  scale_color_manual(values = c("blue","black","red"))+
  theme_clean(base_size = 12)+
  geom_text(aes(MRscore,ORR+0.01, label=Types))+ 
  geom_abline(slope=0.18, intercept=0.27,color="gray",size=1)
ggplot(data1,aes(MRscore,ORR))+
  geom_point(aes(size=N,color=OS))+
  scale_color_manual(values = c("blue","black","red"))+
  theme_stata(base_size = 12)+
  geom_text(aes(MRscore,ORR+0.02, label=Types))+ 
  geom_abline(slope=0.415, intercept=0.249,color="gray",size=1)    
