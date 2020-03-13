##GC_datasets: GSE62254, GSE57303, GSE84437, GSE15459, GSE26253, GSE29272

# my GC IO data
library(data.table)
library(tidyverse)
options(stringsAsFactors = F)
expr<-fread("../dataset/IO_dataset/GC_arranged/GC_expression_tpm.txt",header = T)%>%as.data.frame()
meta<-fread("../dataset/IO_dataset/GC_arranged/GC_clinicaldata.txt",header = T)%>%as.data.frame()
colnames(expr)<-gsub("X","GC",colnames(expr))
colnames(expr)[1]<-"Gene"
rownames(meta)<-paste0("GC",meta$Subject.ID)
colnames(meta)[1]<-"sampleID"
meta$sampleID<-paste0("GC",meta$sampleID)
expr[1:3,1:3]
meta[1:3,1:3]
keypath<-fread("../dataset/TCGA_data/Nod_TLR_genes.csv")[,c(2,4)]
DE<-fread("../dataset/TCGA_data/MRgenes_CancerTYpes_DEGtag.csv")
keypathwaygeneExpr<-merge(keypath,DE,by="Gene")
df<-subset(keypathwaygeneExpr,abs(Log2FC)>=5)
head(df)
keygene<-df$Gene
head(keypath)
GC_keygene<-merge(keypath,expr,by="Gene")
GC_keygene<-GC_keygene[,-2]
GC_keygene<-GC_keygene[!which(duplicated(GC_keygene$Gene)),]
library(pheatmap)
mat<-data.frame(row.names = GC_keygene$Gene,GC_keygene[,-1])
mat1<-mat[which(rownames(mat)%in%df$Gene),]
mat2<-mat[which(rownames(mat)%in%test),]
dfcol<-data.frame(row.names = meta$sampleID,Efficacy=meta$Efficacy)
pheatmap(scale(mat,center = F),border_color = NA,annotation_col = dfcol)
pheatmap(scale(log(mat1+1)),border_color = NA,scale = "column",
         cluster_cols = F,
         annotation_col = dfcol)
pheatmap(scale(log(mat2+1)),border_color = NA,scale = "column",
         cluster_cols = T,
         annotation_col = dfcol)
Heatmap(scale(log(mat1+1)))
library(ComplexHeatmap)
library(survival)
library(survminer)
library(survivalROC)
library(forestplot)
Tmat<-as.data.frame(t(mat))
Tmat1<-as.data.frame(t(mat1))
Tmat$sampleID<-rownames(Tmat)
Tmat1$sampleID<-rownames(Tmat1)
sv_data<-merge(meta,Tmat,by="sampleID")
sv_data1<-merge(meta,Tmat1,by="sampleID")
surv_cut_CCL5 <- surv_cutpoint(
  sv_data,
  time = "PFS",
  event = "Status.1",
  variables = c("CCL5")
)
summary(surv_cut_CCL5)
plot(surv_cut_CCL5)
sv_data$CCL5_group<-sample(c("High","Low"),nrow(sv_data),replace = T)
sv_data$CCL5_group<-ifelse(sv_data$CCL5>=25.53985,'High',"Low")
sv_data$CCL2_group<-sample(c("High","Low"),nrow(sv_data),replace = T)
sv_data$CCL2_group<-ifelse(sv_data$CCL2>=mean(sv_data$CCL2),'High',"Low")
sv_data$DEFA1_group<-sample(c("High","Low"),nrow(sv_data),replace = T)
sv_data$DEFA1_group<-ifelse(sv_data$DEFA1>=mean(sv_data$DEFA1),'High',"Low")
sv_data$CXCL9_group<-sample(c("High","Low"),nrow(sv_data),replace = T)
sv_data$CXCL9_group<-ifelse(sv_data$CXCL9>=mean(sv_data$CXCL9),'High',"Low")
sv_data$CXCL10_group<-sample(c("High","Low"),nrow(sv_data),replace = T)
sv_data$CXCL10_group<-ifelse(sv_data$CXCL10>=mean(sv_data$CXCL10),'High',"Low")
sv_data$DEFA3_group<-sample(c("High","Low"),nrow(sv_data),replace = T)
sv_data$DEFA3_group<-ifelse(sv_data$DEFA3>=mean(sv_data$DEFA3),'High',"Low")
fit1<-survfit(Surv(OS,Status) ~ CCL5_group,
              data = sv_data)
fit1<-survfit(Surv(OS,Status) ~ DEFA3_group,
              data = sv_data)
fit_CCL5_pfs<-survfit(Surv(PFS,Status.1) ~ CCL5_group,
                      data = sv_data)
sv_data$Efficacy1<-sample(c("Response","Non_response"),nrow(sv_data),replace = T)
sv_data$Efficacy1<-ifelse(sv_data$Efficacy=="PR"|sv_data$Efficacy=="SD","Response","Non_response")
fit_CCL5_pfs
p1<-ggsurvplot( fit_CCL5_pfs,
                data=sv_data,
                risk.table = TRUE,
                pval = TRUE,
                title = "PFS",
                palette = c("blue","red"),
                #facet.by = "Efficacy",
                legend.title="CCL5",
                risk.table.col = "strata",
                surv.median.line = "hv",
                risk.table.y.text.col = T,
                risk.table.y.text = FALSE )
p1
#C-index
library(tableone)
## PS matching
library(Matching)
## Weighted analysis
library(survey)
library(reshape2)
vars<-keygene
## Construct a table
tabUnmatched <- CreateTableOne(vars = colnames(sv_data)[-c(1:21)], strata = "CCL5_group", data = sv_data, test = FALSE)
## Show table with SMD
print(tabUnmatched, smd = TRUE)
summary(tabUnmatched, smd = TRUE)
test<-c("CCL5","CARD9","CASP4","CXCL2","HSP90AA1","MFN1","NEK7","RHOA","RNF31")
boothit=200
for(i in 1:boothit){
  case=subset(sv_data,CCL5_group=="Low")
  control=subset(sv_data,CCL5_group=="High")
  bootindex.case=sample(1:nrow(case),replace=T)
  boot.case.data=case[bootindex.case,]
  bootindex.control=sample(1:nrow(control),replace=T)
  boot.control.data=control[bootindex.control,]
  boot.data=rbind(boot.case.data,boot.control.data)
  dstr.boot=svydesign(id=~1, 
                      #prob=~inv_weight, 
                      #fpc=~ssize, 
                      data=boot.data)
  boot.fit=svycoxph(Surv(sv_data$OS,sv_data$Status) ~CCL5, data=boot.data, x=TRUE,design=dstr.boot)
  library(Hmisc)
  library(survival)
  library(survcomp)
  library(survey)
  cindex.train=1-rcorr.cens(lp.boot,Surv(boot.data$time, boot.data$Statuso))
  cindex.test=1-rcorr.cens(lp_=.test,Surv(boot.data$time,boot.data$Statuso))
  bias=rep(1,bootit)
  bias=abs(cindex.train-cindex.test)
} 

#commparisons
my_comparisons<-list(c("PD", "SD"), c("PD", "PR"), c("PR", "SD"))
fun_to_plot <- function(data, group, variable) {
  p <- ggviolin(data, x=group, y=variable,fill = group, 
                palette = c("#00AFBB", "#E7B800", "#FC4E07"), 
                add = "jitter", shape=group)+
    stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means(label.y = 125)
  return(p)
}
fun_to_plot(sv_data,"Efficacy","CCL5")
fun_to_plot(sv_data,"Efficacy","CCL2")
####wilcox 检验批处理示例
#使用其中的 summaryBy() 以方便按分组计算均值、中位数
library(doBy)
expr[1:3,1:3]
which(duplicated(expr$Genesymbol))
expr=expr[-c(13689,13692),]
gene<-data.frame(row.names = expr$Gene,expr[,-1])
group<-data.frame(sample=sv_data$sampleID,group=sv_data$Efficacy1)
head(group)
result <- NULL
for (n in 1:nrow(gene)) {
  gene_n <- data.frame(t(gene[n,]))
  gene_id <- names(gene_n)[1]
  names(gene_n)[1] <- 'gene'
  gene_n$sample <- rownames(gene_n)
  gene_n <- merge(gene_n, group, by = 'sample', all.x = TRUE)
  
  gene_n$group <- factor(gene_n$group)
  p_value <- wilcox.test(gene~group, gene_n)$p.value
  if (!is.na(p_value) & p_value < 0.05) {
    stat <- summaryBy(gene~group, gene_n, FUN = c(mean, median))
    result <- rbind(result, c(gene_id, as.character(stat[1,1]), stat[1,2], stat[1,3], as.character(stat[2,1]), stat[2,2], stat[2,3], p_value))
  }
}
result <- data.frame(result)
names(result) <- c('Gene', 'Non_response', 'mean_NR', 'median_NR', 'Response', 'mean_R', 'median_R', 'p_value')
result$Log2FC<-log(as.numeric(result$mean_R)/as.numeric(result$mean_NR))
write.table(result, '../dataset/IO_dataset/GC_arranged/gene.wilcox.txt', sep = '\t', row.names = FALSE, quote = FALSE)
##MRscore function
expr<-fread("../dataset/IO_dataset/GC_arranged/GC_expression_tpm.txt",header = T)%>%as.data.frame()
meta<-fread("../dataset/IO_dataset/GC_arranged/GC_clinicaldata.txt",header = T)%>%as.data.frame()
colnames(expr)<-gsub("X","GC",colnames(expr))
colnames(expr)[1]<-"Gene"
rownames(meta)<-paste0("GC",meta$Subject.ID)
colnames(meta)[1]<-"sampleID"
meta$sampleID<-paste0("GC",meta$sampleID)
DEexpr<-fread('../dataset/IO_dataset/GC_arranged/gene.wilcox.txt')

expr2MRscore<-function(expr,DEexpr){
  signature<-fread("../dataset/TCGA_data/annotationRow1.csv")
  MRgene_DEGs<-subset(DEexpr,Gene%in%signature$Gene)
  MRgene_DEGs$Regulation<-sample(c("Up","Down"),nrow(MRgene_DEGs),replace = T)
  MRgene_DEGs$Regulation<-ifelse(MRgene_DEGs$Log2FC>=0,"Up","Down")
  DE_Mgene<-MRgene_DEGs[,c("Gene","Regulation")]
  Total_MRscore<-sum(subset(MRgene_DEGs,Log2FC>=0)$Log2FC)+sum(subset(MRgene_DEGs,Log2FC<=0)$Log2FC)
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

#survival plot function
MRscore2Sva<-function(MRscore,meta,survivalTypes){
  survdata<-merge(MRscore,meta,by="sampleID")
  if(survivalTypes=="OS"){
    surv_cut <- surv_cutpoint(
      survdata,
      time = "OS",
      event = "Status",
      variables = c("MRscore")
    )
    summary(surv_cut)
    survdata$Group<-sample(c("High","Low"),nrow(MRscore),replace = T)
    survdata$Group<-ifelse(survdata$MRscore>=as.numeric(summary(surv_cut)[1]),"High","Low")
    fit<-survfit(Surv(OS,Status) ~ Group,
                 data = survdata)
    p<-ggsurvplot( fit,
                   data=survdata,
                   risk.table = TRUE,
                   pval = TRUE,
                   title = "OS",
                   palette = c("blue","red"),
                   #facet.by = "Efficacy",
                   legend.title="CCL5",
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
                   legend.title="CCL5",
                   risk.table.col = "strata",
                   surv.median.line = "hv",
                   risk.table.y.text.col = T,
                   risk.table.y.text = FALSE )
    
  }
  return(p)
}
MRscore2Sva(MRscore = MRscore,meta = meta,survivalTypes = "PFS")

##MRscore2Heatmap
MRscore2Heatmap<-function(expr,MRscoreGroup){
  require(pheatmap)
  dfcol<-data.frame(row.names = group$sample,Group=as.factor(group$group))
  mat<-subset(expr,expr$Gene%in%DE_Mgene$Gene)
  mat<-data.frame(row.names = mat$Genesymbol,mat[,-1])
  pheatmap(scale(mat,center = F),annotation_col = dfcol,
           cluster_cols = F,
           border_color = NA,
           show_colnames = F)
}
##use TCGA STAD DEexprData to compute the GC_MRscore
TCGA_DE<-fread("../dataset/IO_dataset/GC_arranged/TCGA_degenes_STAD.txt")%>%as.data.frame()
head(TCGA_DE)
colnames(TCGA_DE)[c(1,5)]<-c("Gene","Log2FC")
TCGA_DE<-subset(TCGA_DE,abs(Log2FC)>=5)
head(DEexpr)
MRscore_TCGA<-expr2MRscore(expr = expr,DEexpr = TCGA_DE)
MRscore<-expr2MRscore(expr = expr,DEexpr = DEexpr)
MRscore[[1]]
MRscore[[2]]
MRscore_TCGA[[1]]
MRscore_TCGA[[2]]
head(meta)
MRscore_TCGA_meta<-merge(meta,MRscore_TCGA[[2]][,c(1,4)],by="sampleID")
MRscore_TCGA_meta$Efficacy1<-sample(c("Response","Non_response"),nrow(MRscore_TCGA_meta),replace = T)
MRscore_TCGA_meta$Efficacy1<-ifelse(MRscore_TCGA_meta$Efficacy=="PR"|MRscore_TCGA_meta$Efficacy=="SD","Response","Non_response")
MRscore_meta<-merge(meta,MRscore[[2]][,c(1,4)],by="sampleID")
MRscore_meta$Efficacy1<-sample(c("Response","Non_response"),nrow(MRscore_meta),replace = T)
MRscore_meta$Efficacy1<-ifelse(MRscore_meta$Efficacy=="PR"|MRscore_meta$Efficacy=="SD","Response","Non_response")
fun_to_plot(MRscore_TCGA_meta,"Efficacy","MRscore")
fun_to_plot(MRscore_TCGA_meta,"Efficacy1","MRscore")
fun_to_plot(MRscore_meta,"Efficacy","MRscore")
fun_to_plot(MRscore_meta,"Efficacy1","MRscore")
surv_cut <- surv_cutpoint(
  MRscore_meta,
  time = "PFS",
  event = "Status.1",
  variables = c("MRscore")
)
surv_cut <- surv_cutpoint(
  MRscore_meta,
  time = "OS",
  event = "Status",
  variables = c("MRscore")
)
summary(surv_cut)
plot(surv_cut)
MRscore_TCGA_meta$MRscore_group<-sample(c("High","Low"),nrow(sv_data),replace = T)
MRscore_TCGA_meta$MRscore_group<-ifelse(MRscore_TCGA_meta$MRscore>=2.03,'High',"Low")
MRscore_meta$MRscore_group<-sample(c("High","Low"),nrow(sv_data),replace = T)
MRscore_meta$MRscore_group<-ifelse(MRscore_meta$MRscore>=-18.2,'High',"Low")
fit_OS<-survfit(Surv(OS,Status) ~ MRscore_group,
                      data = MRscore_meta)
p1<-ggsurvplot( fit_pfs,
                data=MRscore_meta,
                risk.table = TRUE,
                pval = TRUE,
                title = "OS",
                palette = c("blue","red"),
                #facet.by = "Efficacy",
                legend.title="CCL5",
                risk.table.col = "strata",
                surv.median.line = "hv",
                risk.table.y.text.col = T,
                risk.table.y.text = FALSE )
p1

#DEseq using the counts from RSEM 
library(DESeq2)
library(data.table)
library(tidyverse)
options(stringsAsFactors = F)
gc_counts<-fread("../dataset/IO_dataset/GC_arranged/byRSEM/GC_IO_count.txt")
gc_counts[1:10,1:4]
gc_counts$gene_id<-gsub("_.*","",gc_counts$gene_id)
library(clusterProfiler)
library(org.Hs.eg.db)
gene_ID=select(org.Hs.eg.db,keys = gc_counts$gene_id,column = "SYMBOL",keytype = "ENSEMBL" ,multiVals = "first")
head(gene_ID)
colnames(gc_counts)[1]<-"ENSEMBL"
gc_counts<-merge(gene_ID,gc_counts,by="ENSEMBL")
gc_counts<-gc_counts[,-1]
colnames(gc_counts)[1]<-"Gene"
gc_counts<-gc_counts%>%group_by(Gene)%>%summarise_all(mean)%>%as.data.frame()
gc_counts[1:10,1:4]
gc_counts<-gc_counts[-25787,]
gc_counts[1:10,1:4]
colnames(gc_counts)<-gsub("N.*","",colnames(gc_counts))
write.csv(gc_counts,"../dataset/IO_dataset/GC_arranged/byRSEM/gc_counts_geneSymbol.csv",row.names = F)
meta<-fread("../dataset/IO_dataset/GC_arranged/GC_clinicaldata.txt",header = T)%>%as.data.frame()
colnames(meta)[4]<-"sampleID"
meta$sampleID<-gsub("W.*","",meta$sampleID)
colnames(gc_counts)<-gsub("N.*","",colnames(gc_counts))
gc_counts[1:10,1:4]
meta[1:3,1:4]
index<-which(meta$sampleID%in%colnames(gc_counts))
trans_meta<-meta[index,-c(1:3)]
trans_meta[1:5,1:5]
trans_meta$Group<-trans_meta$Efficacy
trans_meta$Group<-ifelse(trans_meta$Group=="PD"|trans_meta$Group=="SD","R","NR")
gc_counts<-as.data.frame(gc_counts)
gc_counts[1:10,1:4]
gc_counts1<-data.frame(row.names =gc_counts$Gene,gc_counts[,-1])
#gc_counts=as.data.frame(lapply(gc_counts,as.numeric))
gc_counts1[1:10,1:4]
#rownames(GBMcounts)<-rownames(t(GBM))
#colnames(GBMcounts)
dim(gc_counts1)
gc_counts1[1:3,1:3]
GCcol<-data.frame(trans_meta[,c(1,19)])
colnames(GCcol)<-c("ID","Group")
head(GCcol)
Group<-factor(GCcol$Group,levels = c("R","NR"))
str(gc_counts1)
str(GCcol)
str(Group)
write.csv(trans_meta,"../dataset/IO_dataset/GC_arranged/byRSEM/trans_meta.csv",row.names = F)
colData <- data.frame(row.names=colnames(gc_counts1),Group=Group)
head(colData)
gc_counts<-apply(gc_counts1, 2, function(x){round(x,0)})
gc_counts[1:5,1:5]
DEseq_files<-list(countData = gc_counts,
                  colData = colData,
                  design = ~ Group)
save(DEseq_files,file = "../dataset/IO_dataset/GC_arranged/byRSEM/DEseq_files.Rdata")
dds <- DESeqDataSetFromMatrix(countData = gc_counts,
                              colData = colData,
                              design = ~ Group)
dds_filter <- dds[rowSums(counts(dds))>1, ]
nrow(dds_filter)
dds_out <- DESeq(dds_filter)
res <- results(dds_out)
summary(res)
table(res$padj<0.05)
res_deseq <- res[order(res$padj),]
diff_gene_deseq2 <- subset(res_deseq, padj<0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
diff_gene_deseq2 <- row.names(diff_gene_deseq2)
res_diff_data <- merge(as.data.frame(res),as.data.frame(counts(dds_out,normalize=TRUE)),by="row.names",sort=FALSE)
colnames(res_diff_data)[1]<-"Gene"
write.csv(res_diff_data,file = "../dataset/IO_dataset/GC_arranged/byRSEM/DEseq.csv",row.names = F)
DEGs<-fread("../dataset/IO_dataset/GC_arranged/byRSEM/DEseq.csv")[,1:7]
DEGs_filter<-subset(DEGs,abs(log2FoldChange)>=1.5|pvalue<=0.05)
library(EnhancedVolcano)
library(factoextra)
DEGs<-data.frame(row.names = DEGs$Gene,logFC=DEGs$log2FoldChange,
                 AveExpr=DEGs$baseMean,t=DEGs$lfcSE,P.Value=DEGs$pvalue,
                 adj.P.Val=DEGs$padj)
head(DEGs)
VolcanoPlot<-EnhancedVolcano(DEGs,
                             lab = rownames(DEGs),
                             x = "logFC",
                             y = "P.Value",
                             selectLab = rownames(DEGs)[1:5],
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
VolcanoPlot
head(DEGs_filter)
colnames(DEGs_filter)[3]<-"Log2FC"
exprdata<-fread("../dataset/IO_dataset/GC_arranged/byRSEM/DEseq.csv")[,-c(2:7)]
exprdata[1:5,1:5]
expr2MRscore<-function(expr,DEexpr){
  signature<-fread("../dataset/TCGA_data/annotationRow1.csv")
  MRgene_DEGs<-subset(DEexpr,Gene%in%signature$Gene)
  MRgene_DEGs$Regulation<-sample(c("Up","Down"),nrow(MRgene_DEGs),replace = T)
  MRgene_DEGs$Regulation<-ifelse(MRgene_DEGs$Log2FC>=0,"Up","Down")
  DE_Mgene<-MRgene_DEGs[,c("Gene","Regulation")]
  Total_MRscore<-sum(subset(MRgene_DEGs,Log2FC>=0)$Log2FC)+sum(subset(MRgene_DEGs,Log2FC<=0)$Log2FC)
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
GC_MRscore<-expr2MRscore(expr = exprdata,DEexpr = DEGs_filter)
#Total Score = -7.70030276111791
#survival plot function
MRscore<-GC_MRscore[[2]]
head(MRscore)
meta=trans_meta
MRscore2Sva<-function(MRscore,meta,survivalTypes){
  survdata<-merge(MRscore,meta,by="sampleID")
  if(survivalTypes=="OS"){
    surv_cut <- surv_cutpoint(
      survdata,
      time = "OS",
      event = "Status",
      variables = c("MRscore")
    )
    summary(surv_cut)
    survdata$Group<-sample(c("High","Low"),nrow(MRscore),replace = T)
    survdata$Group<-ifelse(survdata$MRscore>=as.numeric(summary(surv_cut)[1]),"High","Low")
    fit<-survfit(Surv(OS,Status) ~ Group,
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
MRscore2Sva(MRscore = MRscore,meta = meta,survivalTypes = "PFS")
MRscore2Sva(MRscore = MRscore,meta = meta,survivalTypes = "OS")
MRscore2Heatmap<-function(expr,DE_Mgene,groupdata){
  require(pheatmap)
  dfcol<-data.frame(row.names = groupdata$sampleID,Group=as.factor(groupdata$Group))
  mat<-subset(expr,expr$Gene%in%DE_Mgene$Gene)
  mat<-data.frame(row.names = mat$Gene,mat[,-1])
  pheatmap(scale(mat,center = T),annotation_col = dfcol,
           cluster_cols = T,
           border_color = NA,
           show_colnames = F)
}
groupdata<-data.frame(sampleID=trans_meta$sampleID,Group=trans_meta$Efficacy)
DE_Mgene<-GC_MRscore[[1]]
head(DE_Mgene)
MRscore2Heatmap(expr = exprdata,DE_Mgene = DE_Mgene,groupdata = groupdata)
