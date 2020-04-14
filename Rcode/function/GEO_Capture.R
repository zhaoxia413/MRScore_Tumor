GEO_Capture<-function(GSE_id){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(c("dplyr","GEOquery","Biobase"))
  require(dplyr)
  require(GEOquery)
  require(Biobase)
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
