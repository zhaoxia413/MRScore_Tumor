library(tidyverse)
library(data.table)
library(GEOquery)
#mc1="GSE13015.Rdata"
mc3="GSE40012.Rdata"
mc5="GSE65682.Rdata"#LUSC tumor:5;normal:43
mc6="GSE21802.Rdata"#LUSC tumor:5;normal:43
mc7="GSE27131.Rdata"#
mc8="GSE28750.Rdata"#LUSC tumor:5;normal:43
mc9="GSE42834.Rdata"#LUSC tumor:5;normal:43
mc10="GSE57065.Rdata"#LUSC tumor:5;normal:43
mc11="GSE68310.Rdata"#LUSC tumor:5;normal:43
mc12="GSE69528.Rdata"#LUSC tumor:5;normal:43
mc13="GSE82050.Rdata"#LUSC tumor:5;normal:43
mc14="GSE111368.Rdata"#LUSC tumor:5;normal:43
datalist<-list(mc3,mc5,mc6,mc7,mc8,mc8,mc10,mc11,mc12,mc13,mc14)
GEOid<-gsub(".Rdata","",datalist)
files<-paste0("../dataset/GEOdatabase/GEO_downloaded/",datalist)
head(files)
downloaded<-list.files("../dataset/GEOdatabase/GEO_downloaded/")
which(datalist%in%downloaded)
load(file = files[1])
group<-data.frame(sampleID=metdata$geo_accession,
                   Group=metdata$`sample type:ch1`,
                   Gender=metdata$`gender:ch1`,
                   Batch=rep(GEOid[1],nrow(metdata)))
group%>%group_by(Group)%>%summarise(n())
exprSet<-expma%>%.[,c(1,which(colnames(.)%in%group1$sampleID))]%>%as.data.frame()
exprSet[1:5,1:5]
write.csv(exprSet,paste0("../dataset/GEOdatabase/GEO_matrix_meta/",GEOid[1],"_expr.csv"),row.names = F)
write.csv(group,paste0("../dataset/GEOdatabase/GEO_matrix_meta/",GEOid[1],"_meta.csv"),row.names = F)
