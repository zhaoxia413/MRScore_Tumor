if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("dplyr","GEOquery","Biobase",
                       "limma","EnhancedVolcano","factoextra"))
library(dplyr)
library(GEOquery)
library(Biobase)
library(limma)
library(EnhancedVolcano)
library(factoextra)
GEOtoExpr<-function(GSE_id){
  gset <- getGEO(GSE_id, GSEMatrix =TRUE, AnnotGPL=TRUE )
  message("GEOdataset download finished !")
  platformID<-levels(factor(gset[[1]]$platform_id))
  gpl <- getGEO(platformID,destdir = getwd()) %>%Table()
  gpl<-data.frame(gpl$ID,gpl[,7])
  colnames(gpl)<-c("probeID",'Symbol')
  gpl$Symbol<-gsub("//.*","",gpl$Symbol)
  class(gset)
  class(gset[[1]])
  dim(pData(gset[[1]]))
  metdata<-pData(gset[[1]])
  print(metdata[1:5,1:5])
  expma<-exprs(gset[[1]])
  dim(expma)
  save(metdata,expma,file = paste0(GSE_id,".Rdata"))
  expma<-data.frame(probeID=rownames(expma),expma)
  annoexpr<-merge(gpl,expma,by="probeID")[,-1]
  annoexpr<-annoexpr%>%group_by(Symbol)%>%summarise_all(mean)
  print(annoexpr[1:5,1:5])
  return(list(annoexpr,metdata))
}

data<-GEOtoExpr(GSE_id = "GSE7904")
meta<-data[[2]]
group<-data.frame(sampleID=meta$geo_accession,Group=meta$characteristics_ch1)

##DEG 

meta<-data[[2]]
group<-data.frame(sampleID=meta$geo_accession,Group=meta$characteristics_ch1)
group%>%group_by(Group)%>%summarise(n())
group<-subset(group,Group!="NO")
group$Group<-ifelse(group$Group=="NB","Normal","Tumor")

datalist<-data

GEOtoDEG<-function(datalist,group){
mat<-t(data.frame(row.names =data[[1]]$Symbol,data[[1]][,-1]))[,-1]
expr<-data.frame(sampleID=rownames(mat),mat)
meta_expr<-merge(group,expr,by="sampleID")
##PCA
rownames(meta_expr)<-meta_expr$sampleID
dat<-as.data.frame(data[[1]][-1,])
dat<-data.frame(row.names = dat$Symbol,dat[,-1])
dat<-dat[,which(colnames(dat)%in%group$sampleID)]
res.pca <- prcomp(t(dat), scale = TRUE)
PCA<-fviz_pca_ind(res.pca,
             label = "none",
             habillage = group$Group,
             addEllipses = TRUE,
             ggtheme = theme_minimal(),
             ellipse.type = "confidence")
print(PCA)
design=model.matrix(~factor(group$Group))
colnames(design)=levels(factor(group$Group))
mycompare<-str_c(colnames(design),collapse = "-")
contrast.matrix<-makeContrasts(mycompare,
                               levels = design)
fit=lmFit(dat,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2) 
DEG<-topTable(fit2, coef=1, n=Inf) %>% na.omit()  ## coef比较分组 n基因数
DEGsig<-subset(DEG,P.Value<=0.05&abs(logFC)>=2)
write.table(DEGsig,"DEGsig.txt",quote = F,sep = "\t")
message("DEGtable downloaded !")
# Volcano plot
VolcanoPlot<-EnhancedVolcano(DEG,
                lab = rownames(DEG),
                x = "logFC",
                y = "P.Value",
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
print(VolcanoPlot)
return(list(DEGsig,PCA,VolcanoPlot))
}

DEG<-GEOtoDEG(datalist = data,group = group)
#example
data<-GEOtoExpr(GSE_id = "GSE2109")
meta<-data[[2]]
group<-data.frame(sampleID=meta$geo_accession,Group=meta$characteristics_ch1)
colnames(group)
group%>%group_by(Group)%>%summarise(n())
group<-subset(group,Group=="metaplasia"|Group=="normal")
DEG<-GEOtoDEG(datalist = data,group = group)
