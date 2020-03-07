# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Fri Mar 6 05:14:39 EST 2020

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE89997", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL17586", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "000000000000000111111111111111"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC"))
write.table(tT, file=stdout(), row.names=F, sep="\t")


################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE89997", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL17586", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
gsms <- "000000000000000111111111111111"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  # set group names

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("STS","LTS")

# set parameters and draw the plot
palette(c("#f4dfdf","#dfeaf4", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE89997", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")
##MRscore analysis
library(data.table)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(ggthemes)
library(EnhancedVolcano)
library(ggrepel)
load(file = "../dataset/IO_dataset/data/GSE89997.Rdata")
DEGres<-fread("../dataset/IO_dataset/data/GEO2R processing_GSE89997.txt")%>%as.data.frame()
expma[1:5,1:5]
metdata[1:5,1:5]
DEGres[1:5,1:7]
DEGres$gene_assignment<-str_extract(DEGres$gene_assignment,"ENST.* ")
DEGres$gene_assignment<-sub("ENST.*? //","",DEGres$gene_assignment)
DEGres$gene_assignment<-sub(" //.*","",DEGres$gene_assignment)
DEGres$gene_assignment<-sub(" ","",DEGres$gene_assignment)
DEGres<-DEGres[-which(is.na(DEGres$gene_assignment)),]
colnames(DEGres)[c(1,7)]<-c("probeID","Symbol")
DEGres<-DEGres[,-1]%>%group_by(Symbol)%>%summarise_all(mean)
DEGres[1:5,1:6]
colnames(DEGres)[1]<-"Gene"
gpl <- getGEO("GPL17586",destdir = getwd()) %>%Table()
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
gpl$Symbol<-str_extract(gpl$Symbol,"ENST.* ")
gpl$Symbol<-sub("ENST.*? //","",gpl$Symbol)
gpl$Symbol<-sub(" //.*","",gpl$Symbol)
gpl$Symbol<-sub(" ","",gpl$Symbol)
head(gpl)
expma<-data.frame(probeID=rownames(expma),expma)
expr<-merge(gpl,expma,by="probeID")
expr<-expr[,-1]
expr[1:5,1:5]
expr<-expr[-which(is.na(expr$Symbol)),]
colnames(expr)[1]<-"Gene"
expr[1:5,1:5]
VolcanoPlot<-EnhancedVolcano(DEGres,
                             lab = rownames(DEGres),
                             x = "logFC",
                             y = "P.Value",
                             selectLab = rownames(DEGres)[1:6],
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
DEGres_filter<-subset(DEGres,P.Value<=0.05)
colnames(DEGres_filter)[1]<-"Gene"
DEGres_filter[1:3,1:6]
expr[1:5,1:5]
expr<-expr
DEexpr<-DEGres_filter
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
  Upexpr_average<-data.frame(sampleID = colnames(Upexpr)[-1],Upscore=apply(Upexpr[,-1], 2, mean))
  Downexpr_average<-data.frame(sampleID = colnames(Downexpr)[-1],Downscore=apply(Downexpr[,-1], 2, mean))
  DE_average<-bind_cols(Upexpr_average,Downexpr_average)[,-3]
  MRscore_value<-DE_average%>%mutate(MRscore=Upscore-Downscore)
  return(list(MRgene_DEGs,MRscore_value))
}
#MRscore<-expr2MRscore(expr = expr,DEexpr = DEGres_filter)
meta<-data.frame(sampleID=metdata$geo_accession,Group=metdata$source_name_ch1)
head(meta)
meta$Group<-ifelse(meta$Group=="PDAC tumor_STS cohort","STS","LTS")
head(MRscore_value)
MRdata<-merge(MRscore_value[,-c(2,3)],meta,by="sampleID")
MRgene_expr<-bind_rows(Upexpr,Downexpr)
write.csv(MRdata,"../dataset/IO_dataset/data/MRscore_GSE89997_sample(-12.3).csv",row.names = F)
write.csv(MRgene_expr,"../dataset/IO_dataset/data/MRgene_exprMatrix_GSE89997.csv",row.names = F)
head(MRdata)
#commparisons
library(ggpubr)
my_comparisons<-list(c("STS", "LTS"))
fun_to_plot <- function(data, group, variable) {
  p <- ggviolin(data, x=group, y=variable,fill = group, 
                #palette = c("#00AFBB", "#E7B800", "#FC4E07"), 
                add = "jitter", shape=group)+
    stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means(label.y = max(data[,2]+1))+
    scale_fill_aaas()
  return(p)
}
fun_to_plot(data = MRdata,group = "Group",variable = "MRscore")
library(pheatmap)
row.names(MRgene_expr)<-MRgene_expr$Gene
MRgene_expr[1:3,1:3]
MRgene_expr<-MRgene_expr[,-1]
dfcol<-data.frame(row.names =meta$sampleID,Group=meta$Group )
pheatmap(scale(log(MRgene_expr+1),center  = F),border_color = NA,
         annotation_col = dfcol,cluster_cols = F)

