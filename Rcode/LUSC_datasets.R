library(data.table)
library(sva)
library(tidyverse)
options(stringsAsFactors = F)
file<-list.files("./raw_meta/")
filepath<-paste0("./raw_meta/",file)
filepath[1]
metalist<-list()
file
filename<-gsub("meta.*","",file)
for (i in seq(length(filepath))) {
  metalist[[i]]<-fread(filepath[i])%>%as.data.frame()
  names(metalist)[i]<-filename[i]
}
levels(factor(meta7$Stage))
meta1<-data.frame(sampleID=metalist[[1]]$geo_accession,
                  Stage=metalist[[1]]$characteristics_ch1.3)
meta1$Stage<-meta1$Stage%>%gsub("Stage:IA","T1",.)%>%gsub("Stage:IB","T1",.)%>%
  gsub("Stage:IIA","T2",.)%>%gsub("Stage:IIB","T2",.)%>%
  gsub("Stage:IIIA","T3",.)%>%gsub("Stage:IIIB","T3",.)%>%  
  gsub("Stage:IV","T4",.)
meta5<-data.frame(sampleID=metalist[[5]]$geo_accession,
                  Stage=metalist[[5]]$characteristics_ch1.3)
meta5$Stage<-meta5$Stage%>%gsub("stage: 1A","T1",.)%>%gsub("stage: 1B","T1",.)%>%gsub("stage: 1","T1",.)%>%
  gsub("stage: 2A","T2",.)%>%gsub("stage: 2B","T2",.)%>%
  gsub("stage: 3A","T3",.)%>%gsub("stage: 3B","T3",.)%>%  
  gsub("stage: 4","T4",.)%>%gsub("stage: n/a","NA",.)

meta6<-data.frame(sampleID=metalist[[6]]$geo_accession,
                  Stage=metalist[[6]]$characteristics_ch1.1)
meta6$Stage<-meta6$Stage%>%gsub("Stage: Ia","T1",.)%>%gsub("Stage: Ib","T1",.)
meta7<-data.frame(sampleID=metalist[[7]]$geo_accession,
                  Stage=metalist[[7]]$characteristics_ch1.5,
                  OS1=metalist[[7]]$characteristics_ch1.13,
                  OS2=metalist[[7]]$characteristics_ch1.14,
                  OS3=metalist[[7]]$characteristics_ch1.15,
                  OS4=metalist[[7]]$characteristics_ch1.16)
meta7$Stage<-meta7$Stage%>%gsub("pathological stage: IA","T1",.)%>%gsub("pathological stage: IB","T1",.)%>%gsub("stage: 1","T1",.)%>%
  gsub("pathological stage: II","T2",.)
meta7<-meta7[-c(227:246),]%>%as.data.frame()
meta7$OS<-sample("os",nrow(meta7),replace = T)
meta7$OS[grep("death:",meta7$OS1)]<-meta7$OS1[grep("death:",meta7$OS1)]
meta7$OS[grep("death:",meta7$OS2)]<-meta7$OS2[grep("death:",meta7$OS2)]
meta7$OS[grep("death:",meta7$OS3)]<-meta7$OS3[grep("death:",meta7$OS3)]
meta7$OS<-gsub("death: ","",meta7$OS)
meta7$OS<-ifelse(meta7$OS=="dead",1,0)
meta7$OStime<-sample("time",nrow(meta7),replace = T)
meta7$OStime[grep("months before relapse/censor:",meta7$OS1)]<-meta7$OS1[grep("months before relapse/censor:",meta7$OS1)]
meta7$OStime[grep("months before relapse/censor:",meta7$OS2)]<-meta7$OS2[grep("months before relapse/censor:",meta7$OS2)]
meta7$OStime[grep("months before relapse/censor:",meta7$OS3)]<-meta7$OS3[grep("months before relapse/censor:",meta7$OS3)]
meta7<-meta7[,-c(3:6)]
meta7$OStime<-meta7$OStime%>%gsub("months before relapse/censor: ","",.)
meta7$OStime<-round(as.numeric(meta7$OStime),0)
levels(factor(meta10$Stage))
meta10<-data.frame(sampleID=metalist[[10]]$geo_accession,
                   Stage=metalist[[10]]$characteristics_ch1.3)
meta10$Stage<-meta10$Stage%>%gsub("stage: 1A","T1",.)%>%gsub("stage: 1B","T1",.)%>%gsub("stage: 1","T1",.)%>%
  gsub("stage: 2A","T2",.)%>%gsub("stage: 2B","T2",.)%>%gsub("stage: 2","T2",.)%>%
  gsub("stage: 3A","T3",.)%>%gsub("stage: 3B","T3",.)%>%  
  gsub("stage: 4","T4",.)%>%gsub("stage: n/a","NA",.)

meta11<-data.frame(sampleID=metalist[[11]]$Accession,
                   Stage=metalist[[11]]$`1:7TH TNM STAGE`,
                   OS=metalist[[11]]$`1:DEATH DUE TO CANCER`,
                   OStime=metalist[[11]]$`1:SURVIVAL AFTER SURGERY (DAYS)`)
levels(factor(meta11$Stage))
meta11$Stage<-meta11$Stage%>%gsub("Ia","T1",.)%>%gsub("Ib","T1",.)%>%
  gsub("IIa","T2",.)%>%gsub("IIb","T2",.)%>%gsub("unknown (I or II)","T2",.)
meta11$Stage[51]<-"T2"
meta11$OStime<-round(as.numeric(meta11$OStime)/365,0)
meta11$OS<-ifelse(meta11$OS=="yes",0,1)
meta12<-data.frame(sampleID=metalist[[12]]$geo_accession,
                   OS=metalist[[12]]$characteristics_ch1.4,
                   OStime=metalist[[12]]$characteristics_ch1.11,
                   Stage=metalist[[12]]$characteristics_ch1.7)
levels(factor(meta12$Stage))
meta12$Stage<-gsub("^.*p","",meta12$Stage)
meta12<-meta12[-c(444:462),]
meta12$OS<-gsub("vital_status: ","",meta12$OS)
meta12$OS<-ifelse(meta12$OS=="Alive",0,1)
meta12$OStime<-gsub("months_to_last_contact_or_death: ","",meta12$OStime)
meta12$OStime<-round(as.numeric(meta12$OStime),0)
meta13<-data.frame(sampleID=metalist[[13]]$geo_accession,
                   OS=metalist[[13]]$characteristics_ch1.12,
                   OStime=metalist[[13]]$characteristics_ch1.10,
                   Stage=metalist[[13]]$characteristics_ch1.4)
levels(factor(meta13$Stage))
meta13$Stage<-meta13$Stage%>%gsub("disease_stage: 1","T1",.)%>%
  gsub("disease_stage: 3","T3",.)%>%
  gsub("disease_stage: --","NA",.)
meta13$OStime<-gsub("^.*: ","",meta13$OStime)
meta13$OStime<-round(as.numeric(meta13$OStime),0)
meta13$OS<-gsub("^.*: ","",meta13$OS)%>%as.numeric()
meta13<-meta13[-c(87:96),]
meta14<-data.frame(sampleID=metalist[[14]]$geo_accession,
                   Stage=metalist[[14]]$V49)
levels(factor(meta14$Stage))
meta14$Stage<-meta14$Stage%>%gsub("IA","T1",.)%>%gsub("IB","T1",.)%>%
  gsub("IIA","T2",.)%>%gsub("IIB","T2",.)%>%
  gsub("IIIA","T3",.)%>%gsub("IIIB","T3",.)%>%  
  gsub("IV","T4",.)%>%gsub("I","",.)%>%gsub("II","",.)
ls()[grep("meta[0-9]",ls())]
names(metalist)
finalmetalist<-list(meta1,meta10,meta11,meta12,meta13,meta14,meta5,meta6,meta7)
names(finalmetalist)<-names(metalist)[c(1,10,11,12,13,14,5,6,7)]
lapply(finalmetalist, colnames)
grep("OS",lapply(finalmetalist, colnames))
needOS<-finalmetalist[-c(3,4,5,9)]
NoneedOS<-finalmetalist[c(3,4,5,9)]
needOS1<-lapply(needOS, function(x){
  x<-x%>%mutate(OS=rep("NA",nrow(x)),OStime=rep("NA",nrow(x)))
})
metamerge1<-bind_rows(needOS1)
metamerge2<-bind_rows(NoneedOS)
metamerge<-rbind(metamerge1,metamerge2)
meta<-fread("./pdata.csv")
meta[1:3,1:3]
colnames(meta)[1]<-"sampleID"
meta[1:3,1:3]
metamerge[1:3,1:3]
LUSC_meta<-merge(meta,metamerge,by="sampleID")
expr<-fread("./mergerall.txt")
expr[1:3,1:3]
colnames(expr)[1]<-"Gene"
expr<-as.data.frame(expr)
LUSC_expr<-data.frame(Gene=expr$Gene,expr[,which(colnames(expr)%in%LUSC_meta$sampleID)])
LUSC_9datasets_expr_meta<-list(LUSC_expr=LUSC_expr,LUSC_meta=LUSC_meta)
save(LUSC_9datasets_expr_meta,file = "./LUSC_9datasets_expr_meta.Rdata")

load(file = "../dataset/dataset_alidation/LUSC/LUSC_9datasets_expr_meta.Rdata")
expr<-LUSC_9datasets_expr_meta$LUSC_expr
meta<-LUSC_9datasets_expr_meta$LUSC_meta
expr[1:3,1:3]
meta[1:3,1:10]
meta<-meta[-which(duplicated(meta$sampleID)),]
expr<-data.frame(row.names = expr$Gene,expr[,-1])
group<-data.frame(row.names = meta$sampleID,datasets=meta$GSE,Group=meta$Tissues)
set.seed(123)
index<-sample(c(1:1011),100)
cbdata<-as.matrix(expr)
cbdata1<-as.matrix(expr[,index])
group1<-group[which(rownames(group)%in%colnames(cbdata1)),]
res.pca <- prcomp(t(cbdata1), scale = TRUE)
library(factoextra)
PCA<-fviz_pca_ind(res.pca,
                  label = "none",
                  habillage = group1$datasets,
                  addEllipses = TRUE,
                  ggtheme = theme_minimal(),
                  ellipse.type = "confidence")
print(PCA)
cbdata[1:3,1:3]
dist_mat1<-dist(t(cbdata1))
clustering1 <- hclust(dist_mat1, method = "complete")
plot(clustering1, labels = group1$datasets)
plot(clustering1, labels =group1$Group)
##sva cross datasets
group%>%group_by(Group,datasets)%>%summarise(n())
#group$bachType<-group$datasets
#for (i in c(1:14)) {
# print(i)
#group$bachType<-gsub(batchs[i],i,group$bachType)
#}
head(group)
mod = model.matrix(~as.factor(datasets), data=group)
#calculating the number of surrogate variables to estimate in a model
n.sv = num.sv(cbdata,mod,method="leek")
n.sv
#2
#sv is far more than the number of datasets
#use the random batchTypes to move batch
group$bachType<-sample(c(1:6),nrow(group),replace = T)
combat_edata = ComBat(dat=cbdata, batch=group$bachType)
#barplot
combat_edata1<-combat_edata[,index]
par(cex = 0.7)
n.sample=ncol(combat_edata1)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
boxplot(combat_edata1,col=cols,main="expression value",las=2)
boxplot(cbdata1,col=cols,main="expression value",las=2)
library(factoextra)
res.pca <- prcomp(t(combat_edata1), scale = TRUE)
PCA<-fviz_pca_ind(res.pca,
                  label = "none",
                  habillage = group1$datasets,
                  addEllipses = TRUE,
                  ggtheme = theme_minimal(),
                  ellipse.type = "confidence")
print(PCA)
##DEGs by limma T vs N
library(limma)
options(stringsAsFactors = F)
library(tidyverse)
library(EnhancedVolcano)
meta<-group
head(meta)
design=model.matrix(~factor(meta$Group)+0)
colnames(design)=levels(factor(meta$Group))
mycompare<-str_c(colnames(design),collapse = "-")
contrast.matrix<-makeContrasts(mycompare,
                               levels = design)
fit=lmFit(combat_edata,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2) 
DEG<-topTable(fit2, adjust="fdr",coef=1, n=Inf) %>% na.omit()  ## coef比较分组 n基因数
# Volcano plot
VolcanoPlot<-EnhancedVolcano(DEG,
                             lab = rownames(DEG),
                             x = "logFC",
                             y = "adj.P.Val",
                             selectLab = rownames(DEG)[1:5],
                             xlab = bquote(~Log[2]~ "fold change"),
                             ylab = bquote(~-Log[10]~italic(P)),
                             pCutoff = 0.05,## pvalue阈值
                             FCcutoff = 1,## FC cutoff
                             xlim = c(-5,5),
                             transcriptPointSize = 1.8,
                             transcriptLabSize = 5.0,
                             colAlpha = 1,
                             legend=c("NS","Log2 FC"," p-value",
                                      " p-value & Log2 FC"),
                             legendPosition = "bottom",
                             legendLabSize = 10,
                             legendIconSize = 3.0)
VolcanoPlot
write.csv(DEG,"../dataset/LUSC_further/DEGs.csv",row.names = T)

head(DEG)
DEG$Gene=rownames(DEG)
DEGs_filter<-subset(DEG,abs(logFC)>=1&P.Value<=0.05)
DEGs_filter$Gene=rownames(DEGs_filter)
expr[1:3,1:3]
expr$Gene<-rownames(expr)
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
head(MRscore)
load(file = "../dataset/LUAD_further/LUSC_9datasets_expr_meta.Rdata")
meta<-LUSC_9datasets_expr_meta$LUSC_meta
head(meta)
meta$OStime<-as.numeric(meta$OStime)
meta$OS<-as.numeric(meta$OS)
#survival plot function
MRscore2Sva<-function(MRscore,meta,survivalTypes){
  library(survival)
  library(survminer)
  library(survivalROC)
  survdata<-merge(MRscore,meta,by="sampleID")
  if(survivalTypes=="OS"){
    surv_cut <- surv_cutpoint(
      survdata,
      time = "OStime",
      event = "OS",
      variables = c("MRscore")
    )
    summary(surv_cut)
    survdata$Group<-sample(c("High","Low"),nrow(MRscore),replace = T)
    survdata$Group<-ifelse(survdata$MRscore>=as.numeric(summary(surv_cut)[1]),"High","Low")
    fit<-survfit(Surv(OStime,OS) ~ Group,
                 data = survdata)
    p<-ggsurvplot( fit,
                   data=survdata,
                   risk.table = TRUE,
                   pval = TRUE,
                   title = "OS",
                   palette = c("blue","red"),
                   #facet.by = "Efficacy",
                   legend.title="MIRscore",
                   risk.table.col = "strata",
                   surv.median.line = "hv",
                   risk.table.y.text.col = T,
                   risk.table.y.text = FALSE )
    p}
  else{
    surv_cut <- surv_cutpoint(
      survdata,
      time = "PFS",
      event = "Status.1",
      variables = c("MRscore")
    )
    summary(surv_cut)
    survdata$Group<-sample(c("High","Low"),nrow(MRscore),replace = T)
    survdata$Group<-ifelse(survdata$MRscore>=as.numeric(summary(surv_cut)[1]),"High","Low")
    fit<-survfit(Surv(PFS,Status.1) ~ Group,
                 data = survdata)
    p<-ggsurvplot( fit,
                   data=survdata,
                   risk.table = TRUE,
                   pval = TRUE,
                   title = "PFS",
                   palette = c("blue","red"),
                   #facet.by = "Efficacy",
                   legend.title="MIRscore",
                   risk.table.col = "strata",
                   surv.median.line = "hv",
                   risk.table.y.text.col = T,
                   risk.table.y.text = FALSE )
    
  }
  return(p)
}
MRscore2Sva(MRscore = MRscore,meta = meta,survivalTypes = "OS")
LUSC_9datasets_expr_meta$LUSC_meta<-meta
  
