source("../../../Git/TuMicroAnalysor/R/requirements.R")
source("../../../Git/TuMicroAnalysor/R/setTaxonClass.R")
source("../../../Git/TuMicroAnalysor/R/taxonClassSplit.R")
source("../../../Git/TuMicroAnalysor/R/reads2abundance.R")
source("../../../Git/TuMicroAnalysor/R/OTUtable_merge.R")
source("../../../Git/TuMicroAnalysor/R/compositionPlotDataMake.R")
source("../../../Git/TuMicroAnalysor/R/compositionPlot.R")
library(survminer)
library(survival)
library(survivalROC)
library(forestplot)
options(stringsAsFactors = F)
meta<-fread("../dataset/TCGA_microreads/data/Metadata-TCGA-Kraken-17625-Samples.csv")%>%as.data.frame()
colnames(meta)[1]<-"sampleID"
micro[1:3,1:3]
meta[1:3,1:3]
meta1<-meta[grep("TCGA",meta$filename),]#6827
colnames(meta)
meta$filename<-gsub("^.*TCGA","TCGA",meta$filename)
meta$filename<-str_sub(meta$filename,1,28)
meta1$filename<-gsub("^.*TCGA","TCGA",meta1$filename)
meta1$filename<-str_sub(meta1$filename,1,28)
head(meta1$filename)
colnames(meta1)[3]<-"tcgaID"
levels(factor(meta$experimental_strategy))
dnameta<-subset(meta1,experimental_strategy=="RNA-Seq")
rnameta<-subset(meta1,experimental_strategy=="WGS")
dnameta<-dnameta[-which(duplicated(dnameta$tcgaID)),]#1933
rnameta<-rnameta[-which(duplicated(rnameta$tcgaID)),]#4118
dnameta$tcgaID<-str_sub(dnameta$tcgaID,1,12)
rnameta$tcgaID<-str_sub(rnameta$tcgaID,1,12)
meta$tcgaID<-str_sub(meta$filename,1,12)
bothmeta<-meta[which(meta$tcgaID%in%intersect(dnameta$tcgaID,rnameta$tcgaID)),]#1795
head(bothmeta$sampleID)
head(dnameta$sampleID)
head(rnameta$sampleID)
OTUtable<-fread("../dataset/TCGA_microreads/userProcessedData/Kraken-TCGA-OTUtable-17625-Samples.csv")%>%as.data.frame()
OTUtable[1:3,1:2]#1993 taxonomy,10230 sampels
OTUtable_rna<-OTUtable[,c(1,which(colnames(OTUtable)%in%rnameta$sampleID))]#2631
OTUtable_dna<-OTUtable[,c(1,which(colnames(OTUtable)%in%dnameta$sampleID))]#1445
OTUtable_both<-OTUtable[,c(1,which(colnames(OTUtable)%in%bothmeta$sampleID))]#1227

OTUtablelist<-list(rna=OTUtable_rna,dna=OTUtable_dna,both=OTUtable_both)
metalist<-list(rna=rnameta,dna=dnameta,both=bothmeta)
OTU_meta_info<-list(OTUtablelist=OTUtablelist,metalist=metalist)
save(OTU_meta_info,file = "../dataset/TCGA_microreads/userProcessedData/OTU_meta_info.Rdata")
path="../dataset/TCGA_microreads/userProcessedData/"
OTUtable2surv<-function(OTUtable,meta,platform,forRange){
  sumreads<-data.frame(sampleID=colnames(OTUtable)[-1],sumreads=apply(OTUtable[,-1], 2, sum))
  meta_sumreads<-merge(sumreads,meta,by="sampleID")
  colnames(meta_sumreads)[13]<-"tcgaID"
  meta_sumreads$OS<-meta_sumreads$vital_status
  meta_sumreads$OS<-ifelse(meta_sumreads$OS=="Alive",0,ifelse(meta_sumreads$OS=="Dead",1,"NA"))
  meta_sumreads$OStime<-round(meta_sumreads$days_to_death/365,0)
  meta_sumreads$OS<-as.numeric(meta_sumreads$OS)
  meta_sumreads$OStime<-as.numeric(meta_sumreads$OStime)
  write.csv(meta_sumreads,paste0(path,platform,"meta_sumreads_TCGA_.csv"),row.names = F)
  survdata<-subset(meta_sumreads,OStime!="NA"&OS!="NA")
  surv_cut_OS <- surv_cutpoint(
    meta_sumreads,
    time = "OStime",
    event = "OS",
    variables = c("sumreads")
  )
  summary(surv_cut_OS)
  plot(surv_cut_OS)
  meta_sumreads$Microgroup<-sample(c("High","Low"),nrow(meta_sumreads),replace = T)
  meta_sumreads$Microgroup<-ifelse(meta_sumreads$sumreads>=5277556,"High","Low")
  fit<-survfit(Surv(OStime,OS) ~ Microgroup,
               data = meta_sumreads)
  fit
  p<-ggsurvplot(fit, data=meta_sumreads,
                risk.table = TRUE,
                pval = TRUE,
                ggtheme = theme_survminer(),
                palette = c("blue","red"),
                #facet.by = "Types",
                legend.title="MicroReads",
                risk.table.col = "strata",
                surv.median.line = "hv",
                risk.table.y.text.col = T,
                risk.table.y.text = FALSE )
  p
  colnames(meta_sumreads)[13]<-"Types"
  meta_sumreads$Types<-gsub("TCGA-","",meta_sumreads$Types)
  
  datalist<-split.data.frame(meta_sumreads,f=meta_sumreads$Types,drop = F)
  surv_cut_OS<-list()
  fit<-list()
  p<-list()
  for (i in forRange) {
    print(i)
    surv_cut_OS[[i]] <- surv_cutpoint(
      datalist[[i]],
      time = "OStime",
      event = "OS",
      variables = c("sumreads")
    )
    summary(surv_cut_OS[[i]])
    plot(surv_cut_OS[[i]])
    datalist[[i]]$MIcroReads<-sample(c("High","Low"),nrow(datalist[[i]]),replace = T)
    datalist[[i]]$MIcroReads<-ifelse(datalist[[i]]$sumreads>=summary(surv_cut_OS[[i]])[,1],"High","Low")
    fit[[i]]<-survfit(Surv(OStime,OS) ~ MIcroReads,
                      data = datalist[[i]])
    p[[i]]<-ggsurvplot(fit[[i]],
                       data=datalist[[i]],
                       risk.table = TRUE,
                       pval = TRUE,
                       ggtheme = theme_survminer(),
                       title = paste0("OS_",names(datalist)[i]),
                       palette = c("blue","red"),
                       #facet.by = "Types",
                       legend.title="MicroReads_",
                       risk.table.col = "strata",
                       surv.median.line = "hv",
                       risk.table.y.text.col = T,
                       risk.table.y.text = FALSE )
    pdf(file = paste0("../dataset/TCGA_microreads/userProcessedData/",names(datalist)[i],platform,"_OS.pdf"),width = 5,height = 6,onefile = FALSE)
    print(p[[i]])
    dev.off()
    png(file = paste0("../dataset/TCGA_microreads/userProcessedData/",names(datalist)[i],platform,"_OS.png"))
    print(p[[i]])
    dev.off()
    res=list(meta_sumreads=meta_sumreads,plot=p)
  }
  return(res)
}
length(levels(factor(rnameta$investigation)))# 23cancertypes
length(levels(factor(dnameta$investigation)))# 8 cancertypes
length(levels(factor(bothmeta$investigation)))# 8 cancertypes
survPlot_rna<-OTUtable2surv(OTUtable = OTUtable_rna,meta = rnameta,
                        platform = "_RNA_",forRange = c(1:23))
survPlot_dna<-OTUtable2surv(OTUtable = OTUtable_dna,meta = dnameta,
                        platform = "_DNA_",forRange = c(1:8))
survPlot_both<-OTUtable2surv(OTUtable = OTUtable_both,meta = bothmeta,
                        platform = "_RNA/DNA_",forRange = c(1:8))
microScore<-fread("../dataset/TCGA_data/MRscore_TCGA_Patients9359.csv")%>%as.data.frame()
colnames(microScore)[1]<-"tcgaID"
colnames(microScore)[1]
microScore$tcgaID[1:2]
microScore$tcgaID<-str_sub(microScore$tcgaID,1,12)
microScore_readsCor<-function(meta_sumreads){
  MIRscore_microReads<-merge(microScore,meta_sumreads[,c(1,2,4)],by="tcgaID")
  MIRscore_microReads$sumreads1<-log10(MIRscore_microReads$sumreads)
  plot_cor_with_label<-function(data,x,y,inputFacet){
    require(ggplot2)
    require(ggthemes)
    require(ggsci)
    cor.coef = TRUE
    corr_eqn <- function(x,y, method='pearson', digits = 2) {
      corr_coef <- round(cor.test(x, y, method=method)$estimate, digits = digits)
      corr_pval <- tryCatch(format(cor.test(x,y, method=method)$p.value, 
                                   scientific=TRUE),
                            error=function(e) NA)
      paste(method, 'r = ', corr_coef, ',', 'p =', round(as.numeric(corr_pval),4))
    }
    p<-data%>%
      ggplot(aes(.[,x],log(.[,y]+1)))+
      geom_point(size = 2,alpha=0.3)+
      theme_bw(base_size = 12)+
      ggtitle(label = paste0(x,"-",y))+
      facet_wrap(~Types,scales = "free",ncol = 5)+
      geom_smooth(method = "lm")
    p
    if (cor.coef) {
      resCor <- data.frame(facets = unique(data[, inputFacet]))
      for(i in seq_along(resCor$facets)) {
        foo <- data[data[,inputFacet] == resCor$facets[i], ]
        resCor$text[i] <- corr_eqn(foo[,x], foo[,y])
      }
      colnames(resCor)[1] <- inputFacet
      resCor$text<-gsub("pearson ","",resCor$text)
      resCor1<-resCor
      resCor2<-resCor
      resCor1$text<-gsub(" , p.*$","",resCor1$text)
      resCor2$text<-gsub("^.*p","p",resCor2$text)
      p <- p + geom_text(data = resCor1, color="red",size=4,
                         aes(x=2,
                             max(log(data[,y]+1)-0.05, na.rm = TRUE),
                             label = text)) + 
        geom_text(data = resCor2, color="red",size=4,
                  aes(x=2,
                      1 * max(log(data[,y]+1), na.rm = TRUE),
                      label = text))+
        theme_igray(base_size = 12)
    }
    return(p) 
  }
  p<-plot_cor_with_label(data = MIRscore_microReads,
                         x = "MRscore",
                         y = "sumreads1",
                         inputFacet = "Types")
  print(p)
  return(p)
}
dna_sumreads<-fread("../dataset/TCGA_microreads/userProcessedData/_DNA_meta_sumreads_TCGA_.csv")%>%as.data.frame()
rna_sumreads<-fread("../dataset/TCGA_microreads/userProcessedData/_RNA_meta_sumreads_TCGA_.csv")%>%as.data.frame()
microScore_readsCor(dna_sumreads)
microScore_readsCor(rna_sumreads)
#compare DNA and RNA have no overlap of samples
OTUtable_dna[1:3,1:3]
OTUtable_dna<-OTUtable_dna%>%mutate(Methods=rep("DNA",nrow(.)))
OTUtable_rna<-OTUtable_rna%>%mutate(Methods=rep("RNA",nrow(.)))
renameMicro<-data.frame(Taxon=OTUtable_dna$Taxon,micro=paste0("OTU",c(1:1993)))
df1<-merge(renameMicro,OTUtable_rna,by="Taxon")[,-1]
df2<-merge(renameMicro,OTUtable_dna,by="Taxon")[,-1]
df_rna<-t(data.frame(row.names = df1$micro,df1[,-1]))%>%as.data.frame()
df_dna<-t(data.frame(row.names = df2$micro,df2[,-1]))%>%as.data.frame()
df_rna<-data.frame(sampleID=rownames(df_rna),df_rna)
df_dna<-data.frame(sampleID=rownames(df_dna),df_dna)
metainfo<-select(meta1,c(sampleID,tcgaID,investigation))
dfr<-merge(metainfo,df_rna,by="sampleID")[,-c(1,3)]
dfd<-merge(metainfo,df_dna,by="sampleID")[,-c(1,3)]
dfr$tcgaID<-str_sub(dfr$tcgaID,1,16)
dfd$tcgaID<-str_sub(dfd$tcgaID,1,16)
index<-intersect(dfr$tcgaID,dfd$tcgaID)##314 overlap
dfr_ovlp<-dfr[which(dfr$tcgaID%in%index),]
dfd_ovlp<-dfd[which(dfd$tcgaID%in%index),]
dfr_ovlp<-dfr_ovlp[-which(duplicated(dfr_ovlp$tcgaID)),]
dfr_ovlp<-dfr_ovlp%>%mutate(Method=rep("RNA",nrow(.)))
dfd_ovlp<-dfd_ovlp%>%mutate(Method=rep("DNA",nrow(.)))
df_both<-bind_rows(dfr_ovlp,dfd_ovlp)
data<-melt(df_both,id.vars = c("tcgaID","Method"),variable.name = "OTUs")
head(data)
d1<-subset(data,Method=="RNA")
colnames(d1)[4]<-"RNA_value"
d2<-subset(data,Method=="DNA")
colnames(d2)[4]<-"DNA_value"
OTUtable_ovlp_dna<-dfd_ovlp[,-which(colnames(dfd_ovlp)=="Method")]
OTUtable_ovlp_rna<-dfr_ovlp[,-which(colnames(dfr_ovlp)=="Method")]
OTUtable_ovlp_rna<-data.frame(t(data.frame(row.names = OTUtable_ovlp_rna$tcgaID,OTUtable_ovlp_rna[,-1])))
OTUtable_ovlp_rna$micro<-rownames(OTUtable_ovlp_rna)
OTUtable_ovlp_dna<-data.frame(t(data.frame(row.names = OTUtable_ovlp_dna$tcgaID,OTUtable_ovlp_dna[,-1])))
OTUtable_ovlp_dna$micro<-rownames(OTUtable_ovlp_dna)
head(renameMicro)
OTUtable_ovlp_dna<-merge(renameMicro,OTUtable_ovlp_dna,by="micro")
OTUtable_ovlp_rna<-merge(renameMicro,OTUtable_ovlp_rna,by="micro")
head(metainfo)
groupboth<-metainfo[,c(2,3)]
groupboth$tcgaID<-str_sub(groupboth$tcgaID,1,16)
groupboth$investigation<-gsub("TCGA-","",groupboth$investigation)
colnames(groupboth)<-c("Samples","Group")
head(groupboth)
groupboth<-groupboth[which(groupboth$Samples%in%index),]
groupboth<-groupboth[-which(duplicated(groupboth$Samples)),]
groupboth$Samples<-gsub("-",".",groupboth$Samples)
OTUtable_ovlp_dna<-OTUtable_ovlp_dna[,-1]
OTUtable_ovlp_rna<-OTUtable_ovlp_rna[,-1]
OTUtable_ovlp_dna[1:3,1:3]
OTUtable_ovlp_rna[1:3,1:3]
###compare the composition of DNA and RNA metheds
dir.create("../dataset/TCGA_microreads/byRNA")
dir.create("../dataset/TCGA_microreads/byDNA")
dir.create("../dataset/TCGA_microreads/byOverlap")
dir.create("../dataset/TCGA_microreads/byOverlap/DNA_overlap")
dir.create("../dataset/TCGA_microreads/byOverlap/RNA_overlap")
setwd("../dataset/TCGA_microreads/byRNA")
OTUtable_dna<-OTUtable_dna[,-which(colnames(OTUtable_dna)=="Methods")]
OTUtable_rna<-OTUtable_rna[,-which(colnames(OTUtable_rna)=="Methods")]
taxon_OTU<-setTaxonClass(taxon_matrix = OTUtable_rna)
OTUsplit<-taxonClassSplit(taxon_matrix = taxon_OTU)
OTUsplit$Kindom<-gsub("\\|.*","",OTUsplit$Kindom)
write.csv(OTUsplit,"Kraken-TCGA-OTUsplitReads-RNA2631-Samples",row.names = F)
mergeOTU<-OTUtablemerge(OTUtable_split = OTUsplit)
head(rnameta)
group<-rnameta[c(1,12)]
group$investigation<-gsub("TCGA-","",group$investigation)
head(group)
colnames(group)<-c("Samples","Group")
OTU_abundance<-reads2abundance(OTUmerge = mergeOTU,meta = group)
sapply(OTU_abundance, function(x){sum(x[,-1])})
sapply(OTU_abundance,nrow)#110  440 1993  267  110

compositionPlotData<-compositionPlotDataMake(abundanceList = OTU_abundance,meta = group,
                                             top_nGenus = 20)
head(compositionPlotData)
plotList<-compositionPlot(plotdata = compositionPlotData,
                          ShowTaxonLevels = "Genus")
setwd("../byDNA/")

taxon_OTU<-setTaxonClass(taxon_matrix = OTUtable_dna)
OTUsplit<-taxonClassSplit(taxon_matrix = taxon_OTU)
OTUsplit$Kindom<-gsub("\\|.*","",OTUsplit$Kindom)
write.csv(OTUsplit,"Kraken-TCGA-OTUsplitReads-DNA1463-Samples",row.names = F)
mergeOTU<-OTUtablemerge(OTUtable_split = OTUsplit)
head(dnameta)
group<-dnameta[c(1,12)]
group$investigation<-gsub("TCGA-","",group$investigation)
head(group)
colnames(group)<-c("Samples","Group")
OTU_abundance<-reads2abundance(OTUmerge = mergeOTU,meta = group)
sapply(OTU_abundance, function(x){sum(x[,-1])})
sapply(OTU_abundance,nrow)#110  440 1993  267  110

compositionPlotData<-compositionPlotDataMake(abundanceList = OTU_abundance,meta = group,
                                             top_nGenus = 20)
head(compositionPlotData)
plotList<-compositionPlot(plotdata = compositionPlotData,
                          ShowTaxonLevels = "Genus")
setwd("../dataset/TCGA_microreads/byOverlap/DNA_overlap")
OTUtable_ovlp_dna[,-1]<-apply(OTUtable_ovlp_dna[,-1], 2, as.numeric)
taxon_OTU<-setTaxonClass(taxon_matrix = OTUtable_ovlp_dna)
OTUsplit<-taxonClassSplit(taxon_matrix = taxon_OTU)
OTUsplit$Kindom<-gsub("\\|.*","",OTUsplit$Kindom)
write.csv(OTUsplit,"Kraken-TCGA-OTUsplitReads-DNAovlp314-Samples",row.names = F)
mergeOTU<-OTUtablemerge(OTUtable_split = OTUsplit)
group<-groupboth
head(group)
OTU_abundance<-reads2abundance(OTUmerge = mergeOTU,meta = group)
sapply(OTU_abundance, function(x){sum(x[,-1])})
sapply(OTU_abundance,nrow)#110  440 1993  267  110
compositionPlotData_dna<-compositionPlotDataMake(abundanceList = OTU_abundance,meta = group,
                                             top_nGenus = 20)
head(compositionPlotData)
plotList<-compositionPlot(plotdata = compositionPlotData_dna,
                          ShowTaxonLevels = "Genus")
setwd("../RNA_overlap/")
OTUtable_ovlp_rna[,-1]<-apply(OTUtable_ovlp_rna[,-1], 2, as.numeric)
taxon_OTU<-setTaxonClass(taxon_matrix = OTUtable_ovlp_dna)
OTUsplit<-taxonClassSplit(taxon_matrix = taxon_OTU)
OTUsplit$Kindom<-gsub("\\|.*","",OTUsplit$Kindom)
write.csv(OTUsplit,"Kraken-TCGA-OTUsplitReads-RNAovlp314-Samples",row.names = F)
mergeOTU<-OTUtablemerge(OTUtable_split = OTUsplit)
group<-groupboth
head(group)
OTU_abundance<-reads2abundance(OTUmerge = mergeOTU,meta = group)
sapply(OTU_abundance, function(x){sum(x[,-1])})
sapply(OTU_abundance,nrow)#110  440 1993  267  110
compositionPlotData_rna<-compositionPlotDataMake(abundanceList = OTU_abundance,meta = group,
                                             top_nGenus = 20)
head(compositionPlotData)
plotList<-compositionPlot(plotdata = compositionPlotData_rna,
                          ShowTaxonLevels = "Genus")
setwd("../../../../MRScore_Tumor/")
