library(data.table)
library(dplyr)
options(stringsAsFactors = F)
path<-"../dataset/IO_dataset/GC_arranged/Allclean_Immunotherapy_data/"
allfiles<-list.files("../dataset/IO_dataset/GC_arranged/Allclean_Immunotherapy_data/")
metafiles<-allfiles[grep("clinical",allfiles)]
metafiles
filenames<-gsub("_.*","",metafiles)
filenames[3:4]<-c("Melanoma_phs000452","Melanoma_PRJEB23709")
metalist<-list()
for (i in seq_along(metafiles)) {
  metalist[[i]]<-fread(paste0(path,metafiles[i]))%>%as.data.frame()
  names(metalist)[i]<-filenames[i]
}
head(metalist$GBM)
metalist$GBM$OStime<-round(metalist$GBM$overall.survival..days./365,0)
metalist$GBM$OS<-metalist$GBM$vital.status
metalist$GBM$OS<-ifelse(metalist$GBM$OS=="Alive",0,1)
levels(factor(metalist$GBM$response))#R,N
head(metalist$Hugo)
metalist$Hugo$OStime<-round(metalist$Hugo$overall.survival..days./365,0)
metalist$Hugo$OS<-metalist$Hugo$vital.status
metalist$Hugo$OS<-ifelse(metalist$Hugo$OS=="Alive",0,1)
levels(factor(metalist$Hugo$response))
metalist$Hugo$response<-ifelse(metalist$Hugo$response=="PD","N","R")
head(metalist$Melanoma_PRJEB23709)
metalist$Melanoma_PRJEB23709$OStime<-round(metalist$Melanoma_PRJEB23709$overall.survival..days./365,0)
metalist$Melanoma_PRJEB23709$OS<-metalist$Melanoma_PRJEB23709$vital.status
metalist$Melanoma_PRJEB23709$OS<-ifelse(metalist$Melanoma_PRJEB23709$OS=="Alive",0,1)
levels(factor(metalist$Melanoma_PRJEB23709$response))
metalist$Melanoma_PRJEB23709$response<-ifelse(metalist$Melanoma_PRJEB23709$response=="PD","N","R")
head(metalist$RCC)
metalist$RCC$OStime<-round(metalist$RCC$os_days/365,0)
metalist$RCC$OS<-metalist$RCC$os_censor
levels(factor(metalist$RCC$response))
metalist$RCC$response<-ifelse(metalist$RCC$response=="no clinical benefit","N","R")
head(metalist$Riaz)
metalist$Riaz$OStime<-round(metalist$Riaz$overall.survival..days./365,0)
metalist$Riaz$OS<-metalist$Riaz$vital.status
metalist$Riaz$OS<-ifelse(metalist$Riaz$OS=="Alive",0,1)
levels(factor(metalist$Riaz$response))
metalist$Riaz$response<-ifelse(metalist$Riaz$response=="PD","N","R")
head(metalist$UC)
metalist$UC$OStime<-round(metalist$UC$os,0)
metalist$UC$OS<-metalist$UC$censOS
levels(factor(metalist$UC$`Best Confirmed Overall Response`))
metalist$UC$response<-metalist$UC$`Best Confirmed Overall Response`
metalist$UC$response<-ifelse(metalist$UC$response=="PD","N","R")

countsfiles<-allfiles[grep("count",allfiles)][-c(4:5,7)]
countsfiles
filenames<-gsub("_.*","",countsfiles)
filenames
countslist<-list()
for (i in seq_along(countsfiles)) {
  print(i)
  countslist[[i]]<-fread(paste0(path,countsfiles[i]))%>%as.data.frame()
  names(countslist)[i]<-filenames[i]
  colnames(countslist[[i]])[1]="ENSEMBL"
  countslist[[i]]$ENSEMBL<-gsub("\\..*","",countslist[[i]]$ENSEMBL)
  countslist[[i]]<-countslist[[i]]%>%group_by(ENSEMBL)%>%summarise_all(mean)
}
library(clusterProfiler)
library(org.Hs.eg.db)
gene_ID<-list()
mergedata<-list()
for (i in seq_along(countslist)) {
  gene_ID[[i]]=select(org.Hs.eg.db,keys = countslist[[i]]$ENSEMBL,column = "SYMBOL",keytype = "ENSEMBL" ,multiVals = "first")
  mergedata[[i]]<-merge(gene_ID[[i]],countslist[[i]],by="ENSEMBL")
  countslist[[i]]<-mergedata[[i]][,-1]%>%group_by(SYMBOL)%>%summarise_all(mean)%>%as.data.frame()
  colnames(countslist[[i]])[1]="Gene"
}
countslist$GBM[1:3,1:3]
countslist$Hugo[1:3,1:3]
countslist$Melanoma[1:3,1:3]
countslist$Riaz[1:3,1:3]
metalist$GBM[1:3,1:3]
metalist$Hugo[1:3,1:3]
metalist$Melanoma_PRJEB23709[1:3,1:3]
metalist$Riaz[1:3,1:3]
countslist<-lapply(countslist, function(x){x=x[-25421,]})
metadata<-list(GBM=metalist$GBM,Hugo=metalist$Hugo,Melanoma=metalist$Melanoma_PRJEB23709,Riaz=metalist$Riaz)
IOdatasets=list(metadata=metadata,countsdata=countslist)
save(IOdatasets,file = paste0(path,"IOdatasets.Rdata"))

###DEseq2
path<-"../dataset/IO_dataset/GC_arranged/Allclean_Immunotherapy_data/"
load(file = paste0(path,"IOdatasets.Rdata"))
countslist<-IOdatasets$countsdata
metadata<-IOdatasets$metadata
meta<-list()
index<-list()
dfcol<-list()
Group<-list()
mat<-list()
colData<-list()
dds<-list()
res<-list()
dds_filter<-list()
dds_out<-list()
res_deseq<-list()
res_diff_data<-list()
diff_gene_deseq2<-list()
library(DESeq2)
library(dplyr)
##all the metadata[[3]] are Response so we set the Treatment as the response column
##with the "N" replace "PRE" meaning the pre-treatment
metadata[[3]]$response<-metadata[[i]]$Treatment
metadata[[3]]$response<-ifelse(metadata[[3]]$response=="PRE","N","R")
for (i in seq_along(countslist)) {
  print(i)
  meta[[i]]<-metadata[[i]]
  index[[i]]<-which(colnames(countslist[[i]])%in%meta[[i]]$Sample_ID)
  mat[[i]]<-countslist[[i]][,index[[i]]]
  rownames(mat[[i]])<-countslist[[i]]$Gene
  #mat[[i]]<-data.frame(row.names = countslist[[i]]$Gene,countslist[[i]][,index[[i]]])
  index[[i]]<-which(meta[[i]]$Sample_ID%in%colnames(countslist[[i]]))
  meta[[i]]<-metadata[[i]][index[[i]],]
  dfcol[[i]]<-data.frame(ID=meta[[i]]$Sample_ID,Group=meta[[i]]$response)
  Group[[i]]<-factor(dfcol[[i]]$Group,levels = c("N","R"))
  colData[[i]] <- data.frame(row.names=colnames(mat[[i]]),Group=Group[[i]])
  mat[[i]]<-as.matrix(round(mat[[i]]),0)
  rownames(mat[[i]])<-countslist[[i]]$Gene
  dds[[i]] <- DESeqDataSetFromMatrix(countData = mat[[i]],
                              colData = colData[[i]],
                              design = ~ Group)
  dds_filter[[i]] <- dds[[i]][rowSums(counts(dds[[i]]))>1, ]
  dds_out[[i]] <- DESeq(dds_filter[[i]])
  res[[i]] <- results(dds_out[[i]])
  summary(res[[i]])
  table(res[[i]]$padj<0.05)
  res_deseq[[i]] <- res[[i]][order(res[[i]]$padj),]
  diff_gene_deseq2[[i]] <- subset(res_deseq[[i]], padj<0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
  diff_gene_deseq2[[i]] <- row.names(diff_gene_deseq2[[i]])
  res_diff_data[[i]] <- merge(as.data.frame(res[[i]]),as.data.frame(counts(dds_out[[i]],normalize=TRUE)),by="row.names",sort=FALSE)
  colnames(res_diff_data[[i]])[1]<-"Gene"
  names(res_diff_data)[i]<-names(metadata)[i]
  write.csv(res_diff_data[[i]],paste0("../dataset/IO_dataset/GC_arranged/Allclean_Immunotherapy_data/",names(metadata)[i],"_DEseq.csv"),row.names = F)
}


library(EnhancedVolcano)
library(factoextra)
DEGs<-list()
DEGs_filter<-list()
VolcanoPlot<-list()
for (i in seq_along(res_diff_data)) {
  DEGs[[i]]<-fread(paste0(path,names(metadata)[i],"_DEseq.csv"))[,1:7]
  names(DEGs[i])<-names(metadata)[i]
  DEGs[[i]]<-data.frame(row.names = DEGs[[i]]$Gene,
                        logFC=DEGs[[i]]$log2FoldChange,
                        AveExpr=DEGs[[i]]$baseMean,
                        t=DEGs[[i]]$lfcSE,P.Value=DEGs[[i]]$pvalue,
                        adj.P.Val=DEGs[[i]]$padj)
  VolcanoPlot[[i]]<-EnhancedVolcano(DEGs[[i]],
                               lab = rownames(DEGs[[i]]),
                               x = "logFC",
                               y = "P.Value",
                               selectLab = rownames(DEGs[[i]])[1:5],
                               xlab = bquote(~Log[2]~ "fold change"),
                               ylab = bquote(~-Log[10]~italic(P)),
                               pCutoff = 0.05,## pvalue阈值
                               FCcutoff = 1,## FC cutoff
                               xlim = c(-5,5),
                               ylim = c(0,5),
                               transcriptPointSize = 1.8,
                               transcriptLabSize = 5.0,
                               colAlpha = 1,
                               legend=c("NS","LogFC"," p-value",
                                        " p-value & LogFC"),
                               legendPosition = "bottom",
                               legendLabSize = 10,
                               legendIconSize = 3.0)
  pdf(paste0(path,names(metadata)[i],"DEGs_vol.pdf"))
  print(VolcanoPlot[[i]])
  dev.off()
  png(paste0(path,names(metadata)[i],"DEGs_vol.png"))
  print(VolcanoPlot[[i]])
  dev.off()
  }
##MRscore
DEGs[[1]][1:3,1:3]
DEGs_filter<-list()
for (i in seq_along(DEGs)) {
  DEGs_filter[[i]]<-subset(DEGs[[i]],abs(logFC)>=1&P.Value<=0.05)
  DEGs_filter[[i]]$Gene<-rownames(DEGs_filter[[i]])
  names(DEGs_filter)[i]<-names(metadata)[i]
}
allfiles<-list.files("../dataset/IO_dataset/GC_arranged/Allclean_Immunotherapy_data/")
exprfiles<-allfiles[grep("fpkm",allfiles)]
exprfiles
filenames<-gsub("_.*","",exprfiles[-c(3,5,6)])
filenames
exprfiles<-exprfiles[-c(3,5,6)]
expr<-list()
for (i in seq_along(exprfiles)) {
  expr[[i]]<-fread(paste0(path,exprfiles[i]))%>%as.data.frame()
  names(expr)[i]<-filenames[i]
  colnames(expr[[i]])[1]<-"Gene"
}
expr[[1]][1:3,1:3]
DEGs_filter[[1]][1:3,1:3]
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
MRscorelist<-list()
for (i in seq_along(DEGs_filter)) {
  MRscorelist[[i]]<-expr2MRscore(expr = expr[[i]],DEexpr = DEGs_filter[[i]])
  names(MRscorelist)[[i]]<-names(expr)[i]
}
##servival
library(survival)
library(survminer)
library(survivalROC)
library(forestplot)
MRscore2Sva<-function(MRscore,meta){
  colnames(MRscore)[1]<-"Sample_ID"
  survdata<-merge(MRscore,meta,by="Sample_ID")
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
    return(p)
}
OSplot<-list()
for (i in seq_along(MRscorelist)) {
  OSplot[[i]]<-MRscore2Sva(MRscore =MRscorelist[[i]][[2]],meta = metadata[[i]])
  names(OSplot)[i]<-names(MRscorelist)[i]
  pdf(paste0(path,names(OSplot)[i],"MIR_OS.pdf"))
  print(OSplot[[i]])
  dev.off()
  png(paste0(path,names(OSplot)[i],"MIR_OS.png"))
  print(OSplot[[i]])
  dev.off()
}
IOdatasets_DEGs_MRscore=list(metadata=metadata,countsdata=countslist,DEGs=res_diff_data,
                     DEGs_filter=DEGs_filter,MRscorelist=MRscorelist)
save(IOdatasets_DEGs_MRscore,file = paste0(path,"IOdatasets_DEGs_MRscore.Rdata"))
##explore the mechanism of metadata$Melanoma for the difference of "Pre" and "EDT"
metadata$Melanoma
colnames(MRscorelist[[3]][[2]])[1]<-"Sample_ID"
data<-merge(MRscorelist[[3]][2],metadata$Melanoma,by="Sample_ID")
data$response<-ifelse(data$response=="N","PRE","EDT")
col31<-c("#303841","#D72323","#377F5B","#375B7F","#F2FCFC","#f0027f",
         "#FAF8DE","#666666","#BDF1F6","#023782","#5e4fa2","#F1C40F",
         "#ff7f00","#cab2d6","#240041","#ffff99","#0E3BF0","#a65628",
         "#f781bf","#808FA6","#2EB872","#F0FFE1","#F33535","#011F4E",
         "#82B269","#D3C13E","#3F9DCD","#014E1F","#AFFFDF","#3D002E",
         "#3A554A")
barplot(rep(1,times=length(col31)),col=col31,border=cm.colors(length(col31)),
        axes=FALSE, main="cm.colors"); box()
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(ggsci)
fun_to_plot <- function(data, group, variable,comparisons) {
  p <- ggboxplot(data, x=group, y=variable,fill = group, 
                 #palette = c("#00AFBB", "#E7B800", "#FC4E07"), 
                 add = "jitter")+
    stat_compare_means(comparisons = comparisons,
                       label.y = c(22, 34,36,38,40))+
    scale_fill_aaas()+
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
          axis.title.x = element_blank())
  return(p)
}
colnames(data)[]
levels(factor(data$Treatment))
my_comparisons<-list(c("EDT","PRE"))
p<-fun_to_plot(data,group = "Treatment",variable = "MRscore",comparisons = my_comparisons)
p+facet_wrap(~treatment,scales = "free")
