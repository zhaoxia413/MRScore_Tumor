library(rfUtilities)
library(caret)
library(data.table)
library(tidyverse)
library(dplyr)
load(file = "../dataset/IO_dataset/GC_arranged/Allclean_Immunotherapy_data/IOdatasets_DEGs_MRscore.Rdata")
save(IOdatasets_DEGs_MRscore,file = "../dataset/IO_dataset/GC_arranged/Allclean_Immunotherapy_data/IOdatasets_DEGs_MRscore.Rdata")
MelanomaMRscore<-IOdatasets_DEGs_MRscore$MRscorelist$Melanoma[[2]]
colnames(MelanomaMRscore)[1]<-'sampleID'
Melanoma_meta<-IOdatasets_DEGs_MRscore$metadata$Melanoma
colnames(Melanoma_meta)[2]<-'sampleID'
data<-merge(MelanomaMRscore,Melanoma_meta,by="sampleID")
data$response<-ifelse(data$response=="PD","N","R")
score_hugo<-IOdatasets_DEGs_MRscore$MRscorelist$Hugo[[2]]
meta_hugo<-IOdatasets_DEGs_MRscore$metadata$Hugo
colnames(meta_hugo)[2]<-"sampleID"
hugo<-merge(score_hugo,meta_hugo,by="sampleID")
score_Riaz<-IOdatasets_DEGs_MRscore$MRscorelist$Riaz[[2]]
meta_Riaz<-IOdatasets_DEGs_MRscore$metadata$Riaz
colnames(meta_Riaz)[2]<-"sampleID"
Riaz<-merge(score_Riaz,meta_Riaz,by="sampleID")

score_GBM<-IOdatasets_DEGs_MRscore$MRscorelist$GBM[[2]]
meta_GBM<-IOdatasets_DEGs_MRscore$metadata$GBM
colnames(meta_GBM)[2]<-"sampleID"
colnames(score_GBM)[1]<-"sampleID"
GBM<-merge(score_GBM,meta_GBM,by="sampleID")

keygene_Hugo<-IOdatasets_DEGs_MRscore$MRscorelist$Hugo[[1]]%>%mutate(Gene=rownames(.))%>%mutate(dataset=rep("Hugo",nrow(.)))
keygene_Melanoma<-IOdatasets_DEGs_MRscore$MRscorelist$Melanoma[[1]]%>%mutate(Gene=rownames(.))%>%mutate(dataset=rep("Melanoma",nrow(.)))
keygene_Riaz<-IOdatasets_DEGs_MRscore$MRscorelist$Riaz[[1]]%>%mutate(Gene=rownames(.))%>%mutate(dataset=rep("Riaz",nrow(.)))
keygene_GBM<-IOdatasets_DEGs_MRscore$MRscorelist$GBM[[1]]%>%mutate(Gene=rownames(.))%>%mutate(dataset=rep("GBM",nrow(.)))
IOkeygene<-bind_rows(list(keygene_Hugo,
                          keygene_Melanoma,
                          keygene_Riaz,
                          keygene_GBM))
write.csv(IOkeygene,"IOkeygene.csv",row.names = F)
mergedata<-bind_rows(list(subset(data,select = c(sampleID,MRscore,response)),
                          subset(hugo,select = c(sampleID,MRscore,response)),
                          subset(Riaz,select = c(sampleID,MRscore,response))
                          
                          ,subset(GBM,select = c(sampleID,MRscore,response))
                          ))
library(randomForest)
library(pROC)
library(data.table)
library(caret)
set.seed(1000)
mergedata$response<-as.factor(mergedata$response)
trainIndex<-sample(nrow(mergedata),nrow(mergedata)*0.8)
trainData<-mergedata[trainIndex,]
testData<-mergedata[-trainIndex,]

fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)
rdaFit <- train(response~MRscore, data = trainData, 
                method = "rda", 
                trControl = fitControl, 
                tuneLength = 4,
                metric = "ROC")


RFfit<-randomForest(response~MRscore, trControl = fitControl,
                    data = trainData,ntree=300)
#pre_ran<-predict(rdaFit,newdata = testData)
pre_ran <- predict(RFfit,newdata=testData)
probValue<-data.frame(sampleID=testData$sampleID,predict(RFfit, type = "prob",
                                                         newdata = testData ))
probValue$Probobility<-probValue$R
probValue$Probobility<-ifelse(probValue$Probobility>=0.5,probValue$Probobility,-probValue$N)
probValue$sampleID<-reorder(probValue$sampleID,probValue$Probobility,decreasing = F)
write.csv(mergedata,"MRscore_response_data.csv",row.names = F)
samindex<-sample(c(1:262),50)
df<-mergedata[samindex,]
write.csv(df[,c(1,2)],"MRscore_response_test.csv",row.names = F)
library(ggthemes)
library(ggsci)
p<-ggplot(probValue,aes(sampleID,Probobility,fill=as.factor(Probobility<=0)))+
  geom_col()+
  theme_wsj(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_discrete(name="Probobility",
                      labels=c("Response", "Non-response"))
p


obs_p_ran = data.frame(prob=pre_ran,obs=testData$response)
#输出混淆矩阵
table(testData$response,pre_ran,dnn=c("True","Predicted"))
roc <- roc(testData$response,as.numeric(pre_ran),
           ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE)
round(auc(roc),2)
plot(roc, print.auc=TRUE, colorize = T,
     auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), 
     max.auc.polygon=TRUE,auc.polygon.col="skyblue", 
     print.thres=TRUE,main='RF_ROC')

plot(roc, lwd = 3, col = "red",print.auc=T, colorize = F,
     auc.polygon=T,
     max.auc.polygon=T,
     auc.polygon.col=alpha("lightcoral",alpha=0.1),
     print.thres=T,main='RF_ROC')

library(ggpubr)
library(ggsci)
library(ggthemes)
fun_to_plot <- function(data, group, variable,comparisons,subtype=NULL) {
  p <- ggboxplot(data, x=group, y=variable,fill = group,facet.by = subtype, 
                 #palette = c("#00AFBB", "#E7B800", "#FC4E07"), 
                 add = "jitter")+
    stat_compare_means(comparisons = comparisons,
                       label.y = c(22, 34,36,38,40))+
    scale_fill_aaas()+
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
          axis.title.x = element_blank())
  return(p)
}

levels(factor(data$Treatment))
levels(factor(data$response))
levels(factor(data$Response))
levels(factor(data$treatment))
my_comparisons<-list(c("EDT","PRE"))
my_comparisons1<-list(c("CR","PD"),
                      c("CR","PR"),
                      c("CR","SD"),
                      c("SD","PR"),
                      c("PR","PD"),
                      c("SD","PD"))
my_comparisons2<-list(c("NR","R"))
my_comparisons3<-list(c("Ipilimumab + Nivolumab","Ipilimumab + Pembrolizumab"),
                      c("Ipilimumab + Nivolumab","Nivolumab"),
                      c("Ipilimumab + Nivolumab","Pembrolizumab"),
                      c("Ipilimumab + Pembrolizumab","Nivolumab"),
                      c("Ipilimumab + Pembrolizumab","Pembrolizumab"),
                      c("Nivolumab" ,"Pembrolizumab"))
data$Treatment<-factor(data$Treatment,levels = c('PRE',"EDT"))

p<-fun_to_plot(data,group = "Treatment",variable = "MRscore",comparisons = my_comparisons)
p
p<-fun_to_plot(data,group = "response",variable = "MRscore",comparisons = my_comparisons1,subtype = "Treatment")
p
p<-fun_to_plot(data,group = "Response",variable = "MRscore",comparisons = my_comparisons2,subtype = "Treatment")
p
p<-fun_to_plot(data,group = "treatment",variable = "MRscore",comparisons = my_comparisons3,subtype = "Treatment")
p
p<-fun_to_plot(data,group = "treatment",variable = "MRscore",comparisons = my_comparisons3)
p

