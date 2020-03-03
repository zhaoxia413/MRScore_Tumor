library(dplyr)
library(GEOquery)
library(Biobase)
# Note that GSEMatrix=TRUE is the default 
GSE_id="GSE3494"
GSE_id = "GSE7904"
GEOtoExpr<-function(GSE_id){
gset <- getGEO(GSE_id, GSEMatrix =TRUE, AnnotGPL=TRUE )
message("GEOdataset download finished !")
#BiocManager::install("hgu133plus2.db")
#suppressMessages(library(hgu133plus2.db))
#keytypes(hgu133plus2.db)
#读取探针名和基因名的对应关系表
#ids=toTable(hgu133plus2SYMBOL)
class(gset)
class(gset[[1]])##ExpressionSet
##提取第一个数据集的phenodata
dim(pData(gset[[1]]))
metdata<-pData(gset[[1]])
print(metdata[1:5,1:5])
colnames(metdata)##phenodata信息很多，但用得上的很少
##提取第一个表达矩阵
expma<-exprs(gset[[1]])
dim(expma)
save(metdata,expma,file = paste0(GSE_id,".Rdata"))
expma<-data.frame(probeID=rownames(expma),expma)
platformID<-levels(factor(gset[[1]]$platform_id))
gpl <- getGEO(platformID,destdir = getwd()) %>%Table()
gpl<-data.frame(gpl$ID,gpl$`Gene Symbol`)
colnames(gpl)<-c("probeID",'Symbol')
gpl$Symbol<-gsub("//.*","",gpl$Symbol)
annoexpr<-merge(gpl,expma,by="probeID")[,-1]
# get the mean value of the duplicated genes.
annoexpr<-annoexpr%>%group_by(Symbol)%>%summarise_all(mean)
print(annoexpr[1:5,1:5])
return(list(annoexpr,metdata))
}
load(file ="GPL6480.soft")
x=getGEO(filename = "GPL6480.soft")
#GSE7904:43 tumor, 7 normal breast and 12 normal organelle
data<-GEOtoExpr(GSE_id = "GSE7904")
meta<-data[[2]]
group<-data.frame(sampleID=meta$geo_accession,Group=meta$characteristics_ch1)
##DEG 
library(limma)
library(factoextra)
mat<-t(data.frame(row.names =data[[1]]$Symbol,data[[1]][,-1]))[,-1]
expr<-data.frame(sampleID=rownames(mat),mat)
mat[1:5,1:5]
expr[1:5,1:5]
meta_expr<-merge(group,expr,by="sampleID")
meta_expr[1:5,1:5]
##PCA
rownames(meta_expr)<-meta_expr$sampleID
res.pca <- prcomp(meta_expr[,-c(1:2)], scale = TRUE)
fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             label = "none",
             habillage = meta_expr$Group,
             # col.ind = dfrow, # 颜色对应group信息
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence")
             #legend.title = "Group",## Legend名称
             #repel = TRUE)
##DEG
#GSE7904:43 tumor, 7 normal breast and 12 normal organelle
group<-subset(group,Group!="NO")
group%>%group_by(Group)%>%summarise(n())
group$Group<-ifelse(group$Group=="NB","Normal","Tumor")
dat<-as.data.frame(data[[1]][-1,])
dat<-data.frame(row.names = dat$Symbol,dat[,-1])
dat<-dat[,which(colnames(dat)%in%group$sampleID)]
dat[1:5,1:5]
res.pca <- prcomp(t(dat), scale = TRUE)
fviz_pca_ind(res.pca,
             label = "none",
             habillage = group$Group,
             addEllipses = TRUE,
             ellipse.type = "confidence")
design=model.matrix(~factor(group$Group))
colnames(design)=levels(factor(group$Group))
head(design)
contrast.matrix<-makeContrasts("Tumor-Normal",
                               levels = design)
contrast.matrix
fit=lmFit(dat,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2) 
DEG<-topTable(fit2, coef=1, n=Inf) %>% na.omit()  ## coef比较分组 n基因数
head(DEG)
# Volcano plot
library(EnhancedVolcano)
EnhancedVolcano(DEG,
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
# GPL="GPL97"##下载平台注释
# gpl <- getGEO(GPL,destdir = getwd()) %>% 
# Table()  ##转换为data.frame格式
#  save(gpl,file = "GPL97_annot.Rdata")
# load(file = "GPL97_annot.Rdata")
# head(gpl)
# colnames(gpl)

##取出注释信息
# probe<-gpl %>% 
# select("ID","Gene Symbol","ENTREZ_GENE_ID")
# head(probe)  
# dim(probe)##22283个
