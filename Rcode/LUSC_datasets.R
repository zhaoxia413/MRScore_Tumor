library(data.table)
library(sva)
ex_b_limma[1:3,1:3]
dim(ex_b_limma)
rm(ex_b_limma)
expr<-fread("./mergerall.txt")
meta<-fread("./pdata.csv")
expr[1:3,1:3]
expr<-data.frame(row.names = expr$Symbol,expr[,-1])
group<-data.frame(row.names = meta$GEO,datasets=meta$GSE,Group=meta$Tissues)
set.seed(123)
index<-sample(c(1:2038),100)
cbdata<-as.matrix(expr)
cbdata1<-as.matrix(expr[,index])
res.pca <- prcomp(t(cbdata1), scale = TRUE)
library(factoextra)
PCA<-fviz_pca_ind(res.pca,
                  label = "none",
                  habillage = group1$datasets,
                  addEllipses = TRUE,
                  ggtheme = theme_minimal(),
                  ellipse.type = "confidence")
print(PCA)
group1<-group[which(rownames(group)%in%colnames(cbdata1)),]
cbdata[1:3,1:3]
dist_mat1<-dist(t(cbdata1))
clustering1 <- hclust(dist_mat1, method = "complete")
plot(clustering1, labels = group1$datasets)
plot(clustering1, labels =group1$Group)
##sva cross datasets
library(dplyr)
group%>%group_by(datasets)%>%summarise(n())


#group$bachType<-group$datasets
#for (i in c(1:14)) {
 # print(i)
  #group$bachType<-gsub(batchs[i],i,group$bachType)
#}
head(group)
mod = model.matrix(~as.factor(datasets), data=group)
#calculating the number of surrogate variables to estimate in a model
n.sv = num.sv(cbdata,mod,method="leek")
#3
n.sv = num.sv(cbdata,mod,method="be")
#57
#sv is far more than the number of datasets
#use the random batchTypes to move batch
group$bachType<-sample(c(1:14),nrow(group),replace = T)
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
  write.csv(DEG,"DEGs.csv",row.names = T)
  head(DEG)
  DEG$Gene=rownames(DEG)
  DEGs_filter<-subset(DEG,abs(logFC)>=1&P.Value<=0.05)
  DEGs_filter$Gene=rownames(DEGs_filter)
  expr[1:3,1:3]
  expr$Gene<-rownames(expr)
  expr2MRscore<-function(expr,DEexpr){
    signature<-fread("../../../dataset/TCGA_data/annotationRow1.csv")
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
    Upexpr_average<-data.frame(sampleID = colnames(Upexpr)[-1],Upscore=apply(Upexpr[,-1], 2, mean))
    Downexpr_average<-data.frame(sampleID = colnames(Downexpr)[-1],Downscore=apply(Downexpr[,-1], 2, mean))
    DE_average<-bind_cols(Upexpr_average,Downexpr_average)[,-3]
    MRscore_value<-DE_average%>%mutate(MRscore=Upscore-Downscore)
    return(list(MRgene_DEGs=MRgene_DEGs,MRscore_value=MRscore_value,Total_MRscore=Total_MRscore))
  }  
MRscorelist<-expr2MRscore(expr = expr,DEexpr = DEGs_filter)  
MRscore<-MRscorelist[[2]]
head(MRscore)
levels(factor(group$datasets))
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
view(metalist[[1]])
meta1<-data.frame(sampleID=metalist[[1]]$geo_accession,
                  Stage=metalist[[1]]$characteristics_ch1.3)
view(metalist[[2]])
view(metalist[[3]])
view(metalist[[4]])
view(metalist[[5]])
meta5<-data.frame(sampleID=metalist[[5]]$geo_accession,
                  Stage=metalist[[5]]$characteristics_ch1.3)
view(metalist[[6]])
meta6<-data.frame(sampleID=metalist[[6]]$geo_accession,
                  Stage=metalist[[6]]$characteristics_ch1.1)
view(metalist[[7]])##OS
meta7<-data.frame(sampleID=metalist[[7]]$geo_accession,
                  Stage=metalist[[7]]$characteristics_ch1.5,
                  OS1=metalist[[7]]$characteristics_ch1.13,
                  OS2=metalist[[7]]$characteristics_ch1.14,
                  OS3=metalist[[7]]$characteristics_ch1.15,
                  OS4=metalist[[7]]$characteristics_ch1.16)
view(metalist[[8]])
view(metalist[[9]])
view(metalist[[10]])
meta10<-data.frame(sampleID=metalist[[10]]$geo_accession,
                   Stage=metalist[[10]]$characteristics_ch1.3)
view(metalist[[11]])##OS
meta11<-data.frame(sampleID=metalist[[11]]$Accession,
                   Stage=metalist[[11]]$`1:7TH TNM STAGE`,
                   OS=metalist[[11]]$`1:DEATH DUE TO CANCER`,
                   OStime=metalist[[11]]$`1:SURVIVAL AFTER SURGERY (DAYS)`)
view(metalist[[12]])##OS
meta12<-data.frame(sampleID=metalist[[12]]$geo_accession,
                   OS=metalist[[12]]$characteristics_ch1.4,
                   OStime=metalist[[12]]$characteristics_ch1.11,
                   Stage=metalist[[12]]$characteristics_ch1.7)
view(metalist[[13]])##OS
meta13<-data.frame(sampleID=metalist[[13]]$geo_accession,
                   OS=metalist[[13]]$characteristics_ch1.12,
                   OStime=metalist[[13]]$characteristics_ch1.10,
                   Stage=metalist[[13]]$characteristics_ch1.4)
view(metalist[[14]])
meta14<-data.frame(sampleID=metalist[[14]]$geo_accession,
                   Stage=metalist[[14]]$V49)
