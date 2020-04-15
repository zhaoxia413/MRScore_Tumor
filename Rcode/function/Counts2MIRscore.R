#' require(data.table)
#' require(tidyverse)
#' path2Function="G:/signature/Git/MRScore_Tumor/Rcode/function/"
# mymeta<-fread(paste0(path2Function,"./testdata/counts2DEGs_meta.csv"))%>%as.data.frame()
# mycountsMat<-fread(paste0(path2Function,"./testdata/counts2DEGs_countsMat.csv"))%>%as.data.frame()
# mycountsMat[1:3,1:3]#keytype choose 3: ENSEMBL
test<-Counts2MRscore(countsMat = mycountsMat,meta = mymeta)
Counts2MRscore<-function(countsMat,meta,sigF=1,sigP=0.05,
                         scoreF=1.5,scoreP=0.01){
  source(paste0(path2Function,"./Counts2DEG.R"))
  source(paste0(path2Function,"./expr2MRscore.R"))
  DEGsdata<-counts2DEGs(countsMat = countsMat,
                        meta = meta,
                        DEGsigP = sigF,DEGsigFC =sigP)
  exprdata<-DEGsdata$DEGsResults[,-c(2:7)]
  exprdata[1:5,1:5]
  DEGfilter<-DEGsdata$DEGs_filter
  MRscore<-expr2MRscore(exprMat =exprdata, DEGmat = DEGfilter,
                        scoreF=1.5,scoreP=0.01)
  return(list(DEG_reslut=exprdata[[1]],MRscore=MRscore[[2]]),MRgenes=MRscore[[1]])
}
