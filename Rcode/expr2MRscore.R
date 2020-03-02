##MRscore function
library(survival)
library(survminer)
library(survivalROC)
library(forestplot)
library(data.table)
expr<-fread("../dataset/IO_dataset/GC_arranged/GC_expression_tpm.txt",header = T)%>%as.data.frame()
meta<-fread("../dataset/IO_dataset/GC_arranged/GC_clinicaldata.txt",header = T)%>%as.data.frame()
colnames(expr)<-gsub("X","GC",colnames(expr))
colnames(expr)[1]<-"Gene"
rownames(meta)<-paste0("GC",meta$Subject.ID)
colnames(meta)[1]<-"sampleID"
meta$sampleID<-paste0("GC",meta$sampleID)
DEexpr<-fread('../dataset/IO_dataset/GC_arranged/gene.wilcox.txt')
library()
expr2MRscore<-function(expr,DEexpr){
  signature<-fread("./data/annotationRow1.csv")
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

MRscore<-expr2MRscore(expr = expr,DEexpr = DEexpr)
MRscore[[1]]
MRscore[[2]]
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
