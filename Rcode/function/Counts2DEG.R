# require(data.table)
# require(tidyverse)
# path2Function="G:/signature/Git/MRScore_Tumor/Rcode/function/"
# mymeta<-fread(paste0(path2Function,"./testdata/counts2DEGs_meta.csv"))%>%as.data.frame()
# mycountsMat<-fread(paste0(path2Function,"./testdata/counts2DEGs_countsMat.csv"))%>%as.data.frame()
# mycountsMat[1:3,1:3]#keytype choose 3: ENSEMBL
# head(mymeta)
# test<-counts2DEGs(countsMat = mycountsMat,meta = mymeta)
counts2DEGs<-function(countsMat,meta,
                       DEGsigP=0.05,
                       DEGsigFC=1){
  require(DESeq2)
  require(data.table)
  require(tidyverse)
  require(EnhancedVolcano)
  require(factoextra)
  require(clusterProfiler)
  require(org.Hs.eg.db)
  options(stringsAsFactors = F)
  keys <- keytypes(org.Hs.eg.db)[c(23,2:6,10,25)]
  message("Choose the geneID type in your countsMat:")
  switch(menu(keys) + 1,
         message("Nothing done\n"), 
         GeneKeyType<-keys[1],
         GeneKeyType<-keys[2],
         GeneKeyType<-keys[3],
         GeneKeyType<-keys[4],
         GeneKeyType<-keys[5],
         GeneKeyType<-keys[6],
         GeneKeyType<-keys[7],
         GeneKeyType<-keys[8]
  )
  gene_ID=select(org.Hs.eg.db,keys = countsMat$gene_id,
                 column = "SYMBOL",keytype = GeneKeyType,
                 multiVals = "first")
  head(gene_ID)
  colnames(countsMat)[1]<-GeneKeyType
  countsMat<-merge(gene_ID,countsMat,by=GeneKeyType)
  countsMat<-countsMat[,-1]
  colnames(countsMat)[1]<-"Gene"
  countsMat<-countsMat%>%group_by(Gene)%>%summarise_all(mean)%>%as.data.frame()
  countsMat[1:10,1:4]
  countsMat<-countsMat[-which(is.na(countsMat$Gene)),]
  countsMat1<-data.frame(row.names =countsMat$Gene,countsMat[,-1])
  #countsMat=as.data.frame(lapply(countsMat,as.numeric)
#rownames(GBMcounts)<-rownames(t(GBM))
#colnames(GBMcounts)
  message(paste0("DataInfo: ",dim(countsMat1)[2], " samples & ",dim(countsMat1)[1]," genes"))
 print(countsMat1[1:3,1:3])
colData <- data.frame(row.names=colnames(countsMat1),Group=meta$Group)
message("metaInfo: ")
print(head(colData))
countsMat<-apply(countsMat1, 2, function(x){round(x,0)})
countsMat[1:5,1:5]
DEseq_files<-list(countData = countsMat,
                  colData = colData,
                  design = ~ Group)
save(DEseq_files,file = "./DEseq_files_pre.Rdata")
message("Finish DEseq datalist make")
message("DEseq_files_pre Filelist stored as Rdata")
dds <- DESeqDataSetFromMatrix(countData = countsMat,
                              colData = colData,
                              design = ~ Group)
dds_filter <- dds[rowSums(counts(dds))>1, ]
nrow(dds_filter)
dds_out <- DESeq(dds_filter)
res <- results(dds_out)
summary(res)
table(res$padj<0.05)
res_deseq <- res[order(res$padj),]
diff_gene_deseq2 <- subset(res_deseq, padj<DEGsigP &abs(log2FoldChange) > DEGsigFC)
diff_gene_deseq2 <- row.names(diff_gene_deseq2)
res_diff_data <- merge(as.data.frame(res),as.data.frame(counts(dds_out,normalize=TRUE)),by="row.names",sort=FALSE)
colnames(res_diff_data)[1]<-"Gene"
write.csv(res_diff_data,file = "DEseq_results.csv",row.names = F)
DEGs_filter<-subset(res_diff_data,abs(log2FoldChange)>=DEGsigFC&pvalue<=DEGsigP)
DEGs<-res_diff_data
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
print(VolcanoPlot)
colnames(DEGs_filter)[3]<-"Log2FC"
message("Finish!")
return(list(DEGsResults=res_diff_data,DEGs_filter=DEGs_filter))
}

