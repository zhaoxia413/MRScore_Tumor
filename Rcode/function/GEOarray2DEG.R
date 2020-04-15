## test
data <-GEOtoExpr(GSE_id = "GSE33479")#LUSC tumor:14;normal:13
meta<-data[[2]]
group<-data.frame(sampleID=meta$geo_accession,Group=meta$characteristics_ch1)
exprMat<- data[[1]]
GEOarray2DEG<-function(exprMat,group,DEGsigP=0.05,DEGsigFC=1){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(c("limma","EnhancedVolcano","factoextra"))
  require(limma)
  require(EnhancedVolcano)
  require(factoextra)
  require(dplyr)
  mat<-t(data.frame(row.names =exprMat$Symbol,exprMat[,-1]))[,-1]
  expr<-data.frame(sampleID=rownames(mat),mat)
  meta_expr<-merge(group,expr,by="sampleID")
  ##PCA
  rownames(meta_expr)<-meta_expr$sampleID
  dat<-as.data.frame(exprMat[-1,])
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
  DEGsig<-subset(DEG,P.Value<=DEGsigP&abs(logFC)>=DEGsigFC)
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