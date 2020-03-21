library(data.table)
library(tidyverse)
library(stringr)
micro<-read.csv("../dataset/TCGA_microreads/data/Kraken-TCGA-Raw-Data-17625-Samples.csv")%>%as.data.frame()
meta<-fread("../dataset/TCGA_microreads/data/Metadata-TCGA-Kraken-17625-Samples.csv")%>%as.data.frame()
colnames(micro)[1]<-"sampleID"
colnames(meta)[1]<-"sampleID"
micro[1:3,1:3]
meta[1:3,1:3]
meta1<-meta[grep("TCGA",meta$filename),]
colnames(meta)
meta$filename<-gsub("^.*TCGA","TCGA",meta$filename)
meta$filename<-str_sub(meta$filename,1,12)
meta1$filename<-gsub("^.*TCGA","TCGA",meta1$filename)
meta1$filename<-str_sub(meta1$filename,1,12)
head(meta1$filename)
colnames(meta1)[3]<-"tcgaID"
levels(factor(meta$experimental_strategy))
dnameta<-subset(meta1,experimental_strategy=="RNA-Seq")
rnameta<-subset(meta1,experimental_strategy=="WGS")
dnameta<-dnameta[-which(duplicated(dnameta$tcgaID)),]
rnameta<-rnameta[-which(duplicated(rnameta$tcgaID)),]
bothmeta<-meta[which(meta$filename%in%intersect(dnameta$tcgaID,rnameta$tcgaID)),]
meta3<-meta[,c(1,3,12)]
head(meta3)
colnames(meta3)[c(2,3)]<-c("tcgaID","Types")
levels(factor(meta3$Types))
meta3$Types<-gsub("^.*-","",meta3$Types)
data<-merge(meta3,micro,by="sampleID")
data[9000:9020,1:3]
write.csv("../dataset/TCGA_microreads/data/TCGAmeta_microReads.csv",row.names = F)
source("../../../Git/TuMicroAnalysor/R/requirements.R")
source("../../../Git/TuMicroAnalysor/R/setTaxonClass.R")
source("../../../Git/TuMicroAnalysor/R/taxonClassSplit.R")
source("../../../Git/TuMicroAnalysor/R/reads2abundance.R")
source("../../../Git/TuMicroAnalysor/R/OTUtable_merge.R")
OTUtable<-fread("../dataset/TCGA_microreads/data/Kraken-TCGA-OTUtable-17625-Samples.csv")%>%as.data.frame()
OTUtable[1:3,1:2]
taxon_OTU<-setTaxonClass(taxon_matrix = OTUtable)
OTUsplit<-taxonClassSplit(taxon_matrix = taxon_OTU)
OTUsplit$Kindom<-gsub("\\|.*","",OTUsplit$Kindom)
write.csv(OTUsplit,"Kraken-TCGA-OTUsplitReads-17625-Samples",row.names = F)
mergeOTU<-OTUtablemerge(OTUtable_split = OTUsplit)
head(meta3)
group<-meta3[c(1,3)]
head(group)
colnames(group)<-c("Samples","Group")
OTU_abundance<-reads2abundance(OTUmerge = mergeOTU,meta = group)
sapply(OTU_abundance, function(x){sum(x[,-1])})
sapply(OTU_abundance,nrow)#110  440 1993  267  110
source("../../../Git/TuMicroAnalysor/R/compositionPlotDataMake.R")
compositionPlotData<-compositionPlotDataMake(abundanceList = OTU_abundance,meta = group,
                                             top_nGenus = 20)
head(compositionPlotData)
compositionPlotData<-compositionPlotData%>%subset(.,Group!="GC_Placebo")
source("../../../Git/TuMicroAnalysor/R/compositionPlot.R")
plotList<-compositionPlot(plotdata = compositionPlotData,ShowTaxonLevels = "Genus")

