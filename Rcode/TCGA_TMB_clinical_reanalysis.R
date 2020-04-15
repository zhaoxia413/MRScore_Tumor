library(data.table)
library(survival)
library(survminer)
library(survey)
library(ggpubr)
library(ggsci)
library(ggthemes)
library(pheatmap)
library(ggpubr)
library(ggplot2)
library(ComplexHeatmap)
library(stringr)
library(data.table)
tmb_micro<-fread("../dataset/TCGA_data/MRscore_TMB_patients7832.csv")%>%as.data.frame()
microScore<-fread("../dataset/TCGA_data/MRscore_TCGA_Patients9359.csv")%>%as.data.frame()
clinic1<-fread("../dataset/TCGA_data/TCGA_Clinical_Variables.txt")%>%as.data.frame()
clinic2<-clinic1[c(162,7,48,50,56,22),]
rownames(clinic2)<-clinic2$V1
clinic2<-data.frame(t(clinic2))[-1,]
clinic2$Patient_ID<-rownames(clinic2)
head(clinic2)
clinic<-fread("../dataset/TCGA_data/TCGA_clinical.csv")%>%as.data.frame()
colnames(clinic)[1]<-"Patient_ID"
colnames(microScore)[1]<-"Patient_ID"
clinic$Patient_ID[1:3]
microScore$Patient_ID[1:3]
clinic2$Patient_ID[1:3]
microScore$Patient_ID<-str_sub(microScore$Patient_ID,1,12)
clinic2$Patient_ID<-str_sub(clinic2$Patient_ID,1,12)
microScore$Patient_ID[1:3]
data<-merge(tmb_micro,clinic,by="Patient_ID")
###tmb os mrscore
data$TMB<-data$`Non-silent per Mb`
data$TMB<-ifelse(data$TMB>=10,"TMB_High","TMB_Low")
data%>%group_by(Types,TMB)%>%summarise(n())
fun_to_plot <- function(data, group, variable,facet,comparisons) {
  p <- ggboxplot(data, x=group, y=variable,fill = group, facet.by = facet,font.label = list(size=6),
                 palette = col16[6:12])+
    stat_compare_means(comparisons = comparisons)+
    geom_jitter(col=rgb(0, 0, 00, 0.2),size=1)
  return(p)
}
my_comparisons<-list(c("TMB_High","TMB_Low"))
fun_to_plot(data,"TMB","MRscore","Types",my_comparisons)
 
fun_to_plot(subset(data,Types=="COAD"),"TMB","MRscore","Types",my_comparisons)

require(pROC)
require(randomForest)
require(data.table)
require(caret)
require(ggplot2)
require(ggthemes)
set.seed(1000)
colnames(data)[5]="TMBvalue"
data$OSgroup<-data$OS
data$OSgroup<-ifelse(data$OSgroup==1,"dead","alive")
data$OSgroup<-as.factor(data$OSgroup)

data_high=subset(data,TMB=="TMB_High")
data_low=subset(data,TMB=="TMB_Low")
rfdata=data_high
rfdata=data_low
trainIndex<-sample(nrow(rfdata),nrow(rfdata)*0.8)
trainData<-rfdata[trainIndex,]
testData<-rfdata[-trainIndex,]
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)
RFfit<-randomForest(OSgroup~MRscore, trControl = fitControl,
                    data = trainData,ntree=300,na.action=na.roughfix)
pre_ran <- predict(RFfit,newdata=testData)
obs_p_ran = data.frame(prob=pre_ran,obs=testData$OSgroup)
#输出混淆矩阵
table(testData$OSgroup,pre_ran,dnn=c("True","Predicted"))
#绘制ROC曲线
roc_res <- roc(testData$OSgroup,as.numeric(pre_ran),
                  ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE)
plot(roc_res, print.auc=TRUE, colorize = T,
     auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), 
     max.auc.polygon=TRUE,auc.polygon.col="skyblue", 
     print.thres=TRUE,main='RF_ROC')

data2<-merge(clinic2,data,by="Patient_ID")
head(microScore)
clinic[1:3,1:3]
data1<-merge(microScore,clinic,by="Patient_ID")
data3<-merge(data1,clinic2,by="Patient_ID")

col16<-c("#b3e2cd","#fdcdac","#cbd5e8","#1b9e77","#d95f02","#7570b3",
         "#e7298a","#66a61e","#e6ab02","#a6761d","#666666",
         "#f4cae4","#e6f5c9","#fff2ae","#f1e2cc","#cccccc")
fun_to_plot <- function(data, group, variable,facet,comparisons) {
  p <- ggboxplot(data, x=group, y=variable,fill = group, facet.by = facet,font.label = list(size=6),
                 palette = col16[6:12])+
    stat_compare_means(comparisons = comparisons)+
    geom_jitter(col=rgb(0, 0, 00, 0.2),size=1)
  return(p)
}
my_comparisons<-list(c())
data4<-data3[-grep("[Not Available]",data3$ajcc_tumor_pathologic_pt),]
data4$Stage<-data4$ajcc_tumor_pathologic_pt
data4$Stage<-gsub("T1[a-c]","T1",data4$Stage)
data4$Stage<-gsub("T2[a-c]","T2",data4$Stage)
data4$Stage<-gsub("T3[a-c]","T3",data4$Stage)
data4$Stage<-gsub("T4[a-d]","T4",data4$Stage)
data4<-data4[-which(is.na(data4$Stage)),]
data4<-data4[-which(data4$Stage%in%c("TX","T0")),]
levels(factor(data3$ajcc_tumor_pathologic_pt))
levels(factor(data3$ajcc_nodes_pathologic_pn))
levels(factor(data3$ajcc_metastasis_pathologic_pm))
my_comparisons<-list(c("T1","T2"),
                     c( "T1","T3"),
                     c( "T1","T4"),
                     c("T2","T3"),
                     c("T2","T4"),
                     c( "T4","T3"))
my_comparisons1<-list(c("Early","Late"))
data4$Stage<-factor(data4$Stage,levels = c("T1","T2","T3","T4"))
data4$Stage_group<-data4$Stage
data4$Stage_group<-ifelse(data4$Stage_group%in%c("T1","T2"),"Early","Late")
write.csv(data4,"TCGA_clinical_new.csv",row.names = F)
p<-fun_to_plot(data4,"Stage","MRscore",facet = "Types",comparisons = my_comparisons)
p+theme_few()
p<-fun_to_plot(data4,"Stage","MRscore",facet = NULL,comparisons = my_comparisons)
p
p<-fun_to_plot(data4,"Stage_group","MRscore",facet = "Types",comparisons = my_comparisons1)
p+theme_few()
p<-fun_to_plot(data4,"Stage_group","MRscore",facet = NULL,comparisons = my_comparisons1)
p
##re plot the survival curve of tcga
data1=subset(data4,data4$Stage_group=="Early")
data1=subset(data4,data4$Stage_group=="Late")
<<<<<<< HEAD

datalist<-split.data.frame(data4,f=data4$Types,drop = F)
=======
datalist<-split.data.frame(data1,f=data1$Types,drop = F)
>>>>>>> 8bb5eb81eae04d21b86a863239c44287b94d2c46
df<-datalist[[9]]
colnames(data4)[c(7,8)]<-c("Age","Gender")
df<-subset(data4,Types=="LUAD")
df<-data4
df$years<-df$age_at_initial_pathologic_diagnosis.y
t1<-surv_cutpoint(
  df,
  time = "OS.time",
  event = "OS",
<<<<<<< HEAD
  variables = c("MRscore")
=======
  variables = c("age_at_initial_pathologic_diagnosis.y")
>>>>>>> 8bb5eb81eae04d21b86a863239c44287b94d2c46
)
summary(t1)
t1<-surv_cutpoint(
  df,
  time = "OS.time",
  event = "OS",
  variables = c("MRscore")
)
summary(t1)
df$MRscore_group<-sample(c("High","Low"),nrow(df),replace = T)
df$MRscore_group<-ifelse(df$MRscore>=summary(t1)[,1],"High","Low")
df$Age_group<-sample(c("Mthan50","Lthan50"),nrow(df),replace = T)
<<<<<<< HEAD
df$Age_group<-ifelse(df$Age>=50,"Mthan50","Lthan50")
coxfit1<-coxph(Surv(OS.time,OS) ~ MRscore_group+Stage_group+Age_group,data = df)
=======
df$Age_group<-ifelse(df$age_at_initial_pathologic_diagnosis.y>=50,"Mthan50","Lthan50")
coxfit1<-coxph(Surv(OS.time,OS) ~ MRscore_group+Stage+years,data = df)
>>>>>>> 8bb5eb81eae04d21b86a863239c44287b94d2c46
coxfit2<-coxph(Surv(OS.time,OS) ~ MRscore_group+strata(Stage_group)+Gender,data = df)
coxfit3<-coxph(Surv(OS.time,OS) ~ MRscore_group+strata(Stage_group)+strata(Gender),data = df)
coxfit<-coxph(Surv(OS.time,OS) ~ MRscore_group,data = df)
summary(coxfit1)
summary(coxfit2)
summary(coxfit3)
summary(coxfit)
temp<-cox.zph(coxfit,transform = "km")
plot(temp)
<<<<<<< HEAD
ggforest(coxfit1,main = "Hazard ratio",
=======
ggforest(coxfit,main = "Hazard ratio",
>>>>>>> 8bb5eb81eae04d21b86a863239c44287b94d2c46
         cpositions = c(0.02, 0.22, 0.4), 
         fontsize = 1,
         refLabel = "reference", noDigits = 2)
library(forestmodel)
print(forest_model(coxfit1,factor_separate_line=F, 
                   format_options = list(colour= "black", 
                                         shape = 12, 
                                         text_size = 4, 
                                         banded = T), 
                   theme = theme_forest()))
<<<<<<< HEAD

mutil_res<-summary(coxfit1)
mutil_res_table<-as.data.frame(cbind(mutil_res$coefficients,mutil_res$conf.int[,-c(1,2)]))
rownames(mutil_res_table)<-gsub("low","",rownames(mutil_res$coefficients))
mutil_res_table<-subset(mutil_res_table,`Pr(>|z|)`<=0.05)
write.csv(mutil_res_table,paste0("../dataset/COX_results/mutil_res_table_age_stage_MRscore.csv"))
multi_pick<-rownames(mutil_res_table)
names(multi_pick)<-names(multi_pick)
fmla <- as.formula(paste0("Surv(OS.time,OS) ~",paste0(multi_pick,collapse = '+')))
cox <- coxph(Surv(OS.time,OS) ~ MRscore_group+Stage_group+Age_group, data = df)#构建多因素生存模型
final_res<-summary(cox)
cox=step(cox,direction = "both")
riskScore=predict(cox,type="risk",newdata=df)
risk=as.vector(ifelse(riskScore>1,"high","low"))
names(risk)<-names(multi_pick)
multiCOX_risk_result<-cbind(id=rownames(cbind(df[,c(1,4,5)],riskScore,risk)),
                                 cbind(df[c(1,4,5)],riskScore,risk))
write.csv(multiCOX_risk_result,paste0("../dataset/COX_results/mutation_model-multiCOX_risk_age_stage_MRscore.csv"),quote = T)
rt=read.csv("../dataset/COX_results/mutation_model-multiCOX_risk_age_stage_MRscore.csv",check.names=F,row.names=1)
roc=survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$riskScore, 
                     predict.time =3, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",round(ro$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()
=======
>>>>>>> 8bb5eb81eae04d21b86a863239c44287b94d2c46
#predict newdata woth coxfit1
newdata<-data.frame(ID=c("A1","A2","A3"),
                    MRscore_group=c("High","Low","Low"),
                    Stage=c("T1","T4","T2"),
                    years=c(24,15,60))
sfit<- survfit(coxfit1,newdata=newdata,id=ID)
sfit$n.event

library(forcats)
library(finalfit)
df=df%>%
  mutate(
    status_os=case_when(
      OS==0~0,
      TRUE~1),
      sex=factor(Gender)%>%
        fct_recode("1"="MALE",
                   "0"="FEMALE")%>%
      ff_label("Sex"),
    stage=factor(Stage)%>%
      fct_recode("1"="T1",
                "2" = "T2",
                 "3"="T3",
                 "4"="T4")%>%
      ff_label("Tumor_stage")
  )
survival_object = df %$% Surv(OS.time,OS)
survival_object = df %$% 
  Surv(OS.time/365, OS)
my_survfit = survfit(survival_object ~ 1, data = df)
my_survfit
summary(my_survfit, times = c(0, 1, 2, 3, 4, 5))
dependent_os = "Surv(OS.time/365, OS)"
explanatory = "MRscore_group"
df %>% 
  surv_plot(dependent_os, explanatory, pval = TRUE)
dependent_os = "Surv(OS.time, OS)"
dependent_stage = "Surv(OS.time, Stage)"
dependent_gender = "Surv(OS.time, Gender)"
explanatory = c( "MRscore_group","Stage_group", "Gender","Age_group")

df %>% 
  finalfit(dependent_os, explanatory)
explanatory= c("MRscore_group", "Gender", "strata(Stage)")
df%>% 
  finalfit(dependent_os, explanatory)
df %>% 
  coxphmulti(dependent_os, explanatory) %>% 
  cox.zph() %>% 
  {zph_result <<- .} %>% 
  plot(var=5)


df %>% 
  hr_plot(dependent_os, explanatory)

survdiff(Surv(OS.time,OS)~ MRscore_group +strata(Stage)+strata(Gender),data=df)
fit<-survfit(Surv(OS.time,OS) ~ MRscore_group,ctype = 1,
             conf.type="log-log",
             data = df)
surv_median(fit)
ggsurvplot(fit = fit,facet.by = "Stage",
           data = df,
           pval.method = T,
           risk.table = TRUE,
           pval = TRUE,
           ggtheme = theme_survminer(),
           palette = c("red","blue"),
           #facet.by = "Types",
           legend.title="MRscore_",
           risk.table.col = "strata",
           surv.median.line = "hv",
           risk.table.y.text.col = T,
           risk.table.y.text = FALSE)
ggsurvplot(fit = fit,facet.by = "Stage_group",
           data = df,
           pval.method = T,
           risk.table = TRUE,
           pval = TRUE,
           ggtheme = theme_survminer(),
           palette = c("red","blue"),
           #facet.by = "Types",
           legend.title="MRscore_",
           risk.table.col = "strata",
           surv.median.line = "hv",
           risk.table.y.text.col = T,
           risk.table.y.text = FALSE)
ggsurvplot(fit = fit,facet.by = "Age_group",
           data = df,
           pval.method = T,
           risk.table = TRUE,
           pval = TRUE,
           ggtheme = theme_survminer(),
           palette = c("red","blue"),
           #facet.by = "Types",
           legend.title="MRscore_",
           risk.table.col = "strata",
           surv.median.line = "hv",
           risk.table.y.text.col = T,
           risk.table.y.text = FALSE)

t1<-surv_cutpoint(
  datalist[[9]],
  time = "OS.time",
  event = "OS",
  variables = c("MRscore")
)
summary(t1)
datalist[[9]]$MRscore_group<-sample(c("High","Low"),nrow(datalist[[9]]),replace = T)
datalist[[9]]$MRscore_group<-ifelse(datalist[[1]]$MRscore>=summary(surv_cut_OS[[i]])[,1],"High","Low")
##校正生存曲线
res <- pairwise_survdiff(Surv(OS.time,OS) ~ MRscore_group,data = datalist[[1]])
res$p.value
symnum(res$p.value, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
       symbols = c("****", "***", "**", "*", "+", " "),
       abbr.colnames = FALSE, na = "")

fit<-survfit(Surv(OS.time,OS) ~ MRscore_group,ctype = 1,conf.type="log-log",
                  data = datalist[[9]])

ggsurvplot(fit = fit,data = datalist[[9]],pval.method = T,
           risk.table = TRUE,
           pval = TRUE,
           ggtheme = theme_survminer(),
           palette = c("red","blue"),
           #facet.by = "Types",
           legend.title="MRscore_",
           risk.table.col = "strata",
           surv.median.line = "hv",
           risk.table.y.text.col = T,
           risk.table.y.text = FALSE )
##
library(Rcmdr)
library(RcmdrMisc)

surv_cut_OS<-list()
fit<-list()
p<-list()
###OS
seq_along(datalist)
try(for (i in c(1:20,22,23,25:28)) {
  print(i)
  surv_cut_OS[[i]] <- surv_cutpoint(
    datalist[[i]],
    time = "OS.time",
    event = "OS",
    variables = c("MRscore")
  )
  summary(surv_cut_OS[[i]])
  plot(surv_cut_OS[[i]])
  datalist[[i]]$MRscore_group<-sample(c("High","Low"),nrow(datalist[[i]]),replace = T)
  datalist[[i]]$MRscore_group<-ifelse(datalist[[i]]$MRscore>=summary(surv_cut_OS[[i]])[,1],"High","Low")
  fit[[i]]<-survfit(Surv(OS.time,OS) ~ MRscore_group,
                    data = datalist[[i]])
  p[[i]]<-ggsurvplot(fit[[i]],
                     data=datalist[[i]],
                     risk.table = TRUE,
                     pval = TRUE,
                     ggtheme = theme_survminer(),
                     title = paste0("OS_",names(datalist)[i]),
                     palette = c("red","blue"),
                     #facet.by = "Types",
                     legend.title="MRscore_",
                     risk.table.col = "strata",
                     surv.median.line = "hv",
                     risk.table.y.text.col = T,
                     risk.table.y.text = FALSE )
  pdf(file = paste0("../dataset/TCGA_results/survical_MRscore/",names(datalist)[i],"_OS_Late.pdf"),width = 5,height = 6,onefile = FALSE)
  print(p[[i]])
  dev.off()
  png(file = paste0("../dataset/TCGA_results/survical_MRscore/",names(datalist)[i],"_OS_Late.png"))
  print(p[[i]])
  dev.off()
},silent = T)

##PFS
for (i in c(1:13,15:20,22,23,25:28)) {
  print(i)
  surv_cut_OS[[i]] <- surv_cutpoint(
    datalist[[i]],
    time = "PFS.time",
    event = "PFS",
    variables = c("MRscore")
  )
  summary(surv_cut_OS[[i]])
  plot(surv_cut_OS[[i]])
  datalist[[i]]$MRscore_group<-sample(c("High","Low"),nrow(datalist[[i]]),replace = T)
  datalist[[i]]$MRscore_group<-ifelse(datalist[[i]]$MRscore>=summary(surv_cut_OS[[i]])[,1],"High","Low")
  fit[[i]]<-survfit(Surv(PFS.time,PFS) ~ MRscore_group,
                    data = datalist[[i]])
  p[[i]]<-ggsurvplot(fit[[i]],
                     data=datalist[[i]],
                     risk.table = TRUE,
                     pval = TRUE,
                     ggtheme = theme_survminer(),
                     title = paste0("PFS_",names(datalist)[i]),
                     palette = c("red","blue"),
                     #facet.by = "Types",
                     legend.title="MRscore_",
                     risk.table.col = "strata",
                     surv.median.line = "hv",
                     risk.table.y.text.col = T,
                     risk.table.y.text = FALSE )
  pdf(file = paste0("../dataset/TCGA_results/survical_MRscore/",names(datalist)[i],"_PFS_Late.pdf"),width = 5,height = 6,onefile = FALSE)
  print(p[[i]])
  dev.off()
  png(file = paste0("../dataset/TCGA_results/survical_MRscore/",names(datalist)[i],"_PFS_Late.png"))
  print(p[[i]])
  dev.off()
}
#TMB frequence
library(RColorBrewer)
library(ggridges)

colnames(data)[5]<-"TMB"
subdata<-subset(data,TMB<=8)#data 7827samples subdata 7090
subdata$Types<-subdata$Types[order(subdata$Types,decreasing = F)]
Freq<-ggplot(subdata, 
             aes(TMB,Types, fill = ..density..)) + 
  geom_density_ridges_gradient(aes(height = ..density..),scale = 1,size = 0.3)+#sacle���÷�???
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11,'BrBG')))(32))+
  labs(x="TMB")+
  theme_few(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5))
#facet_wrap(~Immune_Subtype,scales = "free_x")
Freq
Freq1<-ggplot(subdata, 
             aes(MRscore,Types, fill = ..density..)) + 
  geom_density_ridges_gradient(aes(height = ..density..),scale = 1,size = 0.3)+#sacle���÷�???
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(9,'YlGnBu')))(32))+
  labs(x="MRscore")+
  theme_few(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5))
#facet_wrap(~Immune_Subtype,scales = "free_x")
Freq1
##TMB_MRscore_OS_cutoff
datalist1<-split.data.frame(subdata,f=subdata$Types,drop = F)
fun_to_plot <- function(data, group, variable) {
  p <- ggboxplot(data, x=group, y=variable,fill = group, 
                 palette = c("#00AFBB", "#E7B800", "#FC4E07"), 
                 add = "jitter", shape=group)+
    stat_compare_means(comparisons = my_comparisons)
  return(p)
}
surv_cut_OS<-list()
p<-list()
###OS
seq_along(datalist1)
my_comparisons<-list(c("High","Low"))
for (i in seq_along(datalist1)) {
  print(i)
  surv_cut_OS[[i]] <- surv_cutpoint(
    datalist1[[i]],
    time = "OS.time",
    event = "OS",
    variables = c("MRscore")
  )
  summary(surv_cut_OS[[i]])
  plot(surv_cut_OS[[i]])
  datalist1[[i]]$MRscore_level<-sample(c("High","Low"),nrow(datalist1[[i]]),replace = T)
  datalist1[[i]]$MRscore_level<-ifelse(datalist1[[i]]$MRscore>=summary(surv_cut_OS[[i]])[,1],"High","Low")
  datalist1[[i]]$MRscore_level<-factor( datalist1[[i]]$MRscore_level,levels = c("High","Low"))
  p[[i]]<-fun_to_plot(datalist1[[i]],"MRscore_level","TMB")
  pdf(file = paste0("../dataset/TCGA_results/survical_MRscore/",names(datalist1)[i],"_OScutMRscore_TMB.pdf"),width = 5,height = 6,onefile = FALSE)
  print(p[[i]])
  dev.off()
  png(file = paste0("../dataset/TCGA_results/survical_MRscore/",names(datalist1)[i],"_OScutMRscore_TMB.png"))
  print(p[[i]])
  dev.off()
  }
##PFS_cutoff MRscore level
for (i in c(1:11,13:22)) {
  print(i)
  surv_cut_OS[[i]] <- surv_cutpoint(
    datalist1[[i]],
    time = "PFS.time",
    event = "PFS",
    variables = c("MRscore")
  )
  summary(surv_cut_OS[[i]])
  plot(surv_cut_OS[[i]])
  datalist1[[i]]$MRscore_level<-sample(c("High","Low"),nrow(datalist1[[i]]),replace = T)
  datalist1[[i]]$MRscore_level<-ifelse(datalist1[[i]]$MRscore>=summary(surv_cut_OS[[i]])[,1],"High","Low")
  datalist1[[i]]$MRscore_level<-factor( datalist1[[i]]$MRscore_level,levels = c("High","Low"))
  p[[i]]<-fun_to_plot(datalist1[[i]],"MRscore_level","TMB")
  pdf(file = paste0("../dataset/TCGA_results/survical_MRscore/",names(datalist1)[i],"_PFScutMRscore_TMB.pdf"),width = 5,height = 6,onefile = FALSE)
  print(p[[i]])
  dev.off()
  png(file = paste0("../dataset/TCGA_results/survical_MRscore/",names(datalist1)[i],"_PFScutMRscore_TMB.png"))
  print(p[[i]])
  dev.off()
}
plotdata<-bind_rows(datalist1)
p <-ggplot(plotdata,aes(Types,TMB,fill = MRscore_level)) +
  geom_boxplot()+
  scale_fill_manual(values = c("#00AFBB", "#E7B800"))+
  theme_bw()+
  coord_flip()
p  
p <-ggplot(plotdata,aes(Types,TMB,fill = MRscore_level)) +
  geom_boxplot()+
  scale_fill_d3()+
  theme_bw()+
  coord_flip()
p  



