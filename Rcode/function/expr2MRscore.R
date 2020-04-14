#' require(data.table)
#' require(tidyverse)
#' path2Function="G:/signature/Git/MRScore_Tumor/Rcode/function/"
#' source(paste0(path2Function,"./Counts2DEG.R"))
#' mymeta<-fread(paste0(path2Function,"./testdata/counts2DEGs_meta.csv"))%>%as.data.frame()
#' mycountsMat<-fread(paste0(path2Function,"./testdata/counts2DEGs_countsMat.csv"))%>%as.data.frame()
#' DEGsdata<-counts2DEGs(countsMat = mycountsMat,meta = mymeta)
#' exprdata<-DEGsdata$DEGsResults[,-c(2:7)]
#' exprdata[1:5,1:5]
#' DEGfilter<-DEGsdata$DEGs_filter
#' DEGfilter[1:5,1:5]
#' MRscore<-expr2MRscore(exprMat =exprdata, DEGmat = DEGfilter)
expr2MRscore<-function(exprMat,DEGmat,scoreF=1.5,scoreP=0.01){
  require(data.table)
  require(dplyr)
  DEGmat<-subset(DEGmat,abs(Log2FC)>=scoreF&pvalue<=scoreP)
  signature<-fread(paste0(path2Function,"./testdata/signature.csv"))
  MRgene_DEGs<-subset(DEGmat,Gene%in%signature$Gene)
  MRgene_DEGs$Regulation<-sample(c("Up","Down"),nrow(MRgene_DEGs),replace = T)
  MRgene_DEGs$Regulation<-ifelse(MRgene_DEGs$Gene>=0,"Up","Down")
  DE_Mgene<-MRgene_DEGs[,c("Gene","Regulation")]
  Total_MRscore<-sum(subset(MRgene_DEGs,Log2FC>=0)$Log2FC)+sum(subset(MRgene_DEGs,Log2FC<=0)$Log2FC)
  message(paste0("Total Score = ",Total_MRscore))
  upgene<-subset(DE_Mgene,Regulation=="Up")
  downgene<-subset(DE_Mgene,Regulation=="Down")
  Upexpr<-exprMat[which(exprMat$Gene%in%upgene$Gene),]
  Downexpr<-exprMat[which(exprMat$Gene%in%downgene$Gene),]
   if(dim(Upexpr)[1]==0){
     Upexpr_average<-data.frame(sampleID = colnames(Upexpr)[-1],Upscore=rep(0,dim(exprMat)[2]-1))
   } else {Upexpr_average<- data.frame(sampleID = colnames(Upexpr)[-1],Upscore=apply(Upexpr[,-1], 2, mean))}
  if(dim(Downexpr)[1]==0){
    Downexpr_average<-data.frame(sampleID = colnames(Downexpr)[-1],Downscore=rep(0,dim(exprMat)[2]-1))
  }else {
    Downexpr_average<-data.frame(sampleID = colnames(Downexpr)[-1],Downscore=apply(Downexpr[,-1], 2, mean))
  }
  DE_average<-bind_cols(Upexpr_average,Downexpr_average)[,-3]
  MRscore_value<-DE_average%>%mutate(MRscore=Upscore-Downscore)
  message("Finish MRscore compute!")
  return(list(MRgene_FoldChnange=MRgene_DEGs,MRscore_value=MRscore_value))
}

