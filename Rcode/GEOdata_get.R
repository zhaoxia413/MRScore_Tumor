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
library(frma)

library(GEOmetadb)
# ## 查看数据库信息
file.info('GEOmetadb.sqlite')
## 建立连接
con <- dbConnect(SQLite(), 'GEOmetadb.sqlite')
library(dplyr)
db <- src_sqlite('../dataset/GEOdatabase/GEOmetadb.sqlite')
gse_dp <- tbl(db, 'gse')
filter(gse_dp, gse=='GSE2553')
#断开数据库/删除数据库文件
dbDisconnect(con)
file.remove('GEOmetadb.sqlite')
#download rawdata
getGEOSuppFiles("GSE89997")
setwd("GSE89997/")
untar("GSE89997_RAW.tar") ##解压文件
files <- dir(pattern="gz$") ##加载文件
sapply(files, gunzip) ##合并文件
filelist <- dir(pattern="CEL$")
##统一处理CEL文件
library(affy)
library(annotate)
data <- ReadAffy(filenames=filelist)
affydb<-annPkgName(data@annotation,type="db")
require(affydb, character.only=TRUE)
eset<-rma(data,verbose=FALSE)
eset.e <- exprs(eset)  ##得到的矩阵文件
library(annaffy)
symbols<-as.character(aafSymbol(as.character(rownames(eset)),affydb))
genes<-as.character(aafUniGene(as.character(rownames(eset)),affydb))

##download expression matrix
GEOtoExpr<-function(GSE_id){
  gset <- getGEO(GSE_id, GSEMatrix =TRUE, AnnotGPL=TRUE )
  message("GEOdataset download finished !")
  platformID<-levels(factor(gset[[1]]$platform_id))
  gpl <- getGEO(platformID,destdir = getwd()) %>%Table()
  print(head(gpl))
  metaCol<-colnames(gpl)
  message("Choose the gene annotation column:")
    switch(menu(metaCol) + 1,
           message("Nothing done\n"), 
           chooseCol<-metaCol[1],
           chooseCol<-metaCol[2],
           chooseCol<-metaCol[3],
           chooseCol<-metaCol[4],
           chooseCol<-metaCol[5],
           chooseCol<-metaCol[6],
           chooseCol<-metaCol[7],
           chooseCol<-metaCol[8],
           chooseCol<-metaCol[9],
           chooseCol<-metaCol[10],
           chooseCol<-metaCol[11],
           chooseCol<-metaCol[12],
           chooseCol<-metaCol[13],
           chooseCol<-metaCol[14],
           chooseCol<-metaCol[15],
           chooseCol<-metaCol[16]
          )
  gpl<-data.frame(gpl$ID,gpl[,chooseCol])
  colnames(gpl)<-c("probeID",'Symbol')
  gpl$Symbol<-gsub("//.*","",gpl$Symbol)
  gpl<-data.frame(gpl$ID,gpl[,chooseCol])
  colnames(gpl)<-c("probeID",'Symbol')
  #gpl1$Symbol<-str_extract(gpl1$Symbol,"ENST.* ")
  #gpl1$Symbol<-sub("ENST.*? //","",gpl1$Symbol)
  #gpl1$Symbol<-sub(" //.*","",gpl1$Symbol)
  class(gset)
  class(gset[[1]])
  dim(pData(gset[[1]]))
  metdata<-pData(gset[[1]])
  print(metdata[1:5,1:5])
  expma<-exprs(gset[[1]])
  expma<-data.frame(Gene=gpl$probeID,expma)
  dim(expma)
  save(metdata,expma,file = paste0(GSE_id,".Rdata"))
  expma<-data.frame(probeID=rownames(expma),expma)
  annoexpr<-merge(gpl,expma,by="probeID")[,-1]
  annoexpr<-annoexpr%>%group_by(Symbol)%>%summarise_all(mean)
  write.csv(annoexpr,paste0("../dataset/dataset_alidation/",GSE_id,"_expr.csv"),row.names = F)
  write.csv(metdata,paste0("../dataset/dataset_alidation/",GSE_id,"_meta.csv"),row.names = F)
  print(annoexpr[1:5,1:5])
  return(list(annoexpr,metdata))
}
data<-GEOtoExpr(GSE_id = "GSE7904")
##DEG 
meta<-data[[2]]
group<-data.frame(sampleID=meta$geo_accession,Group=meta$characteristics_ch1)
group%>%group_by(Group)%>%summarise(n())
group<-subset(group,Group!="NO")
group$Group<-ifelse(group$Group=="NB","Normal","Tumor")
datalist<-data
GEOdatalist2DEG<-function(datalist,group){
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
GEOdatalist2DEG(datalist = data,group = group)
#example
GSE_id = "GSE3494"#BRCA
df<-GEOtoExpr(GSE_id = "GSE2990")#  BRCA
data<-GEOtoExpr(GSE_id = "GSE2109")#many cancers
data1<-GEOtoExpr(GSE_id = "GSE33479")#LUSC tumor:14;normal:13
data2<-GEOtoExpr(GSE_id = "GSE74706")#LUSC tumor:18;normal:18
data3<-GEOtoExpr(GSE_id = "GSE21933")#LUSC tumor:10;normal:10
data4<-GEOtoExpr(GSE_id = "GSE8569")#LUSC tumor:69;normal:6
data5<-GEOtoExpr(GSE_id = "GSE40275")#LUSC tumor:5;normal:43
data6<-GEOtoExpr(GSE_id ="GSE22863")#LUSC tumor:5;normal:43
data7<-GEOtoExpr(GSE_id ="GSE89997")#PDAC 2 cohort  
##QC
exprSet<-data1
checkQC<-function(exprSet){
  if("arrayQualityMetrics" %in% rownames(installed.packages()) == FALSE) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("arrayQualityMetrics")}
  suppressMessages(library(arrayQualityMetrics))
  arrayQualityMetrics(exprSet, outdir="microarray")  #质控
  browseURL(file.path("microarray", "index.html"))
  MAplot(myData, pairs=TRUE, plot.method="smoothScatter") #MAplot图
  Density<-plotDensity.AffyBatch(exprSet) #密度图
  boxplot<-boxplot(exprSet)               #箱型图
  rnaDeg <- AffyRNAdeg(exprSet)  #查看RNA降解
  plotAffy<-plotAffyRNAdeg(rnaDeg)
  summaryAffyRNAdeg<-summaryAffyRNAdeg(rnaDeg) #获取RNA降解情况
  return(list(Density=Density,
              boxplot=boxplot,
              rnaDeg=rnaDeg,
              plotAffy=plotAffy))
}


meta1<-data1[[2]]
meta2<-data2[[2]]
meta3<-data3[[2]]
meta4<-data4[[2]]
meta5<-data5[[2]]
group1<-data.frame(sampleID=meta1$geo_accession,
                   Group=meta1$`histology:ch1`,
                   Batch=rep("GSE33479",nrow(meta1)))
group1%>%group_by(Group)%>%summarise(n())
group1<-subset(group1,Group=="squamous cell carcinoma"|Group=="normal")
group1$Group<-ifelse(group1$Group=="normal","Normal","LUSC")

group2<-data.frame(sampleID=meta2$geo_accession,
                   Group=meta2$`tissue:ch1`,
                   Batch=rep("GSE74706",nrow(meta2)))
group2%>%group_by(Group)%>%summarise(n())
group2<-subset(group2,Group=="NSCLC"|Group=="Tumor-free lung")
group2$Group<-ifelse(group2$Group=="Tumor-free lung","Normal","LUSC")

group3<-data.frame(sampleID=meta3$geo_accession,
                   Group=meta3$source_name_ch1,
                   Batch=rep("GSE21933",nrow(meta3)))
group3%>%group_by(Group)%>%summarise(n())
group3<-subset(group3,Group=="primary lung tumor"|Group=="primary normal lung tissue")
group3$Group<-ifelse(group3$Group=="primary normal lung tissue","Normal","LUSC")

group4<-data.frame(sampleID=meta4$geo_accession,
                   Group=meta4$characteristics_ch1.1,
                   Batch=rep("GSE8569",nrow(meta4)))
group4%>%group_by(Group)%>%summarise(n())
group4$Group<-ifelse(group4$Group=="Normal lung","Normal","LUSC")

group5<-data.frame(sampleID=meta5$geo_accession,
                   Group=meta5$source_name_ch1,
                   Batch=rep("GSE40275",nrow(meta5)))
group5%>%group_by(Group)%>%summarise(n())
group5<-subset(group5,Group=="normal lung"|Group=="Carcinoma of lung, squamous cell")
group5$Group<-ifelse(group5$Group=="normal lung","Normal","LUSC")

head(group1)
head(group2)
head(group3)
head(group4)
head(group5)
group<-bind_rows(list(group1,group2,group5))
group%>%group_by(Group)%>%summarise(n())
group%>%group_by(Group,Batch)%>%summarise(n())

expr1<-data1[[1]]%>%.[,c(1,which(colnames(.)%in%group1$sampleID))]%>%.[-1,]
expr2<-data2[[1]]%>%.[,c(1,which(colnames(.)%in%group2$sampleID))]%>%.[-1,]
expr3<-data3[[1]]%>%.[,c(1,which(colnames(.)%in%group3$sampleID))]%>%.[-1,]
expr4<-data4[[1]]%>%.[,c(1,which(colnames(.)%in%group4$sampleID))]%>%.[-1,]
expr5<-data5[[1]]%>%.[,c(1,which(colnames(.)%in%group5$sampleID))]%>%.[-1,]
expr1[1:5,1:5]
expr2[1:5,1:5]#log
expr3[1:5,1:5]#2015 genes
expr4[1:5,1:5]#7133 genes
expr5[1:5,1:5]

##

  logNormal<-function(x){
    log(x+1)
  }
library(clusterProfiler)
library(org.Hs.eg.db)
expr5_ID=select(org.Hs.eg.db,keys = expr5$Symbol,column = "SYMBOL",keytype = "ENTREZID" ,multiVals = "first")
colnames(expr5)[1]<-"ENTREZID"
expr5<-merge(expr5_ID,expr5,by="ENTREZID")[,-1]
colnames(expr5)[1]<-"Symbol"
expr5[1:5,1:5]
rm(expr5_ID)
expr1_log<-data.frame(Symbol=expr1[,1],apply(expr1[,-1], 2, logNormal))
expr5_log<-data.frame(Symbol=expr5$Symbol,apply(expr5[,-1], 2, logNormal))
expr<-inner_join(inner_join(expr1,expr2),expr5)
expr<-data.frame(row.names = expr$Symbol,expr[,-1])
expr<-expr[,which(colnames(expr)%in%group$sampleID)]
expr[1:5,1:5]
library(sva)
group%>%group_by(Group,Batch)%>%summarise(n())
group$bachType<-sample(c(1,2,3),nrow(group),replace = T)
group$bachType<-ifelse(group$Batch=="GSE33479",1,ifelse(group$Batch=="GSE40275",2,3))
head(group)
group<-data.frame(row.names = group$sampleID,group)
cbdata<-as.matrix(expr)
dist_mat <- dist(t(cbdata))
clustering <- hclust(dist_mat, method = "complete")
plot(clustering, labels = group$bachType)
plot(clustering, labels =group$Group)
mod = model.matrix(~as.factor(Group), data=group)
n.sv = num.sv(cbdata,mod,method="leek")
combat_edata = ComBat(dat=cbdata, batch=group$bachType, mod=mod, 
                      par.prior=T, prior.plots=T)

par(cex = 0.7)
n.sample=ncol(expr)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
boxplot(expr,col=cols,main="expression value",las=2)
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")
write.csv(expr,"../dataset/dataset_alidation/expr_3LUSC-dataset.csv",row.names = F)
write.csv(exprMat,"../dataset/dataset_alidation/combat_3LUSC-dataset.csv",row.names = F)
write.csv(group,"../dataset/dataset_alidation/combat_3LUSC-dataset_meta.csv",row.names = F)

exprMat<-as.data.frame(combat_edata)
exprMat[1:5,1:5]
head(group)

GEO2DEG<-function(exprMat,group){
  mat<-t(exprMat)
  expr<-data.frame(sampleID=rownames(mat),mat)
  meta_expr<-merge(group[,c(1,2)],expr,by="sampleID")
  ##PCA
  rownames(meta_expr)<-meta_expr$sampleID
  dat<-data.frame(row.names = expr$sampleID,expr[,-1])
  dat<-t(dat)
  dat<-dat[,which(colnames(dat)%in%group$sampleID)]
  res.pca <- prcomp(t(dat), scale = TRUE)
  PCA<-fviz_pca_ind(res.pca,
                    #label = "none",
                    habillage = group$Group,
                    addEllipses = TRUE,
                    ggtheme = theme_minimal(),
                    ellipse.type = "confidence")
  print(PCA)
  ##remove the outlier point: GSM990252,GSM990251,GSM990250,
  ##GSM990249,GSM990248
  outlier<-c("GSM990252","GSM990251","GSM990250","GSM990248","GSM990249")
  dat<-dat[,-which(colnames(dat)%in%outlier)]
  group<-group[-which(group$sampleID%in%outlier),]
  dat<-dat[,which(colnames(dat)%in%group$sampleID)]
  res.pca <- prcomp(t(dat), scale = TRUE)
  PCA<-fviz_pca_ind(res.pca,
                    label = "none",
                    habillage = group$Group,
                    addEllipses = TRUE,
                    ggtheme = theme_minimal(),
                    ellipse.type = "confidence")
  print(PCA)
  
  
  design=model.matrix(~factor(group$Group)+0)
  colnames(design)=levels(factor(group$Group))
  mycompare<-str_c(colnames(design),collapse = "-")
  contrast.matrix<-makeContrasts(mycompare,
                                 levels = design)
  fit=lmFit(dat,design)
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2) 
  DEG<-topTable(fit2, coef=1, n=Inf) %>% na.omit()  ## coef比较分组 n基因数
  colnames(DEG)
  DEGsig<-subset(DEG,adj.P.Val<=0.01&abs(logFC)>=2)
  write.table(DEGsig,"DEGsig.txt",quote = F,sep = "\t")
  message("DEGtable downloaded !")
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
  print(VolcanoPlot)
    return(list(DEGsig,PCA,VolcanoPlot))
}
mc1<-GEOtoExpr(GSE_id = "GSE13015")#LUSC tumor:14;normal:13
mc2<-GEOtoExpr(GSE_id = "GSE20346")#LUSC tumor:18;normal:18
mc3<-GEOtoExpr(GSE_id = "GSE40012")#LUSC tumor:10;normal:10
mc4<-GEOtoExpr(GSE_id = "GSE60244")#LUSC tumor:69;normal:6
mc5<-GEOtoExpr(GSE_id = "GSE65682")#LUSC tumor:5;normal:43
mc6<-GEOtoExpr(GSE_id ="GSE21802")#LUSC tumor:5;normal:43
mc7<-GEOtoExpr(GSE_id ="GSE27131")#
mc8<-GEOtoExpr(GSE_id ="GSE28750")#LUSC tumor:5;normal:43
mc9<-GEOtoExpr(GSE_id ="GSE42834")#LUSC tumor:5;normal:43
mc10<-GEOtoExpr(GSE_id ="GSE57065")#LUSC tumor:5;normal:43
mc11<-GEOtoExpr(GSE_id ="GSE68310")#LUSC tumor:5;normal:43
mc12<-GEOtoExpr(GSE_id ="GSE69528")#LUSC tumor:5;normal:43
mc13<-GEOtoExpr(GSE_id ="GSE82050")#LUSC tumor:5;normal:43
mc14<-GEOtoExpr(GSE_id ="GSE111368")#LUSC tumor:5;normal:43
IO1<-GEOtoExpr(GSE_id ="GSE100797")#Lauss et al., Nat Commun 2017
IO2<-GEOtoExpr(GSE_id ="GSE78220")# Hugo et al., Cell 2016  
IO3<-GEOtoExpr(GSE_id ="GSE93157")#Prat et al., Cancer Res 2017
IO4<-GEOtoExpr(GSE_id = "GSE91061")#Riaz et al., Cell 2017 
library(genefu)
library(breastCancerMAINZ)
library(breastCancerTRANSBIG)
library(breastCancerUPP)
library(breastCancerUNT)
library(breastCancerNKI)
data(breastCancerData)
data.all <- c("transbig7g"=transbig7g, "unt7g"=unt7g, "upp7g"=upp7g,
              "mainz7g"=mainz7g, "nki7g"=nki7g)
dn <- c("transbig", "unt", "upp", "mainz", "nki")
dn.platform <- c("affy", "affy", "affy", "affy", "agilent")
cinfo <- pData(transbig7g)

meta<-data1[[2]]
group<-data.frame(sampleID=meta$geo_accession,Group=meta$characteristics_ch1)
colnames(group)
group%>%group_by(Group)%>%summarise(n())
group<-subset(group,Group=="metaplasia"|Group=="normal")
DEG<-GEOtoDEG(datalist = data,group = group)
