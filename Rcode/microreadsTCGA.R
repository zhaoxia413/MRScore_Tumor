<<<<<<< HEAD
library(data.table)
library(tidyverse)
library(stringr)
micro<-read.csv("../dataset/TCGA_microreads/data/Kraken-TCGA-Raw-Data-17625-Samples.csv")%>%as.data.frame()
meta<-fread("../dataset/TCGA_microreads/data/Metadata-TCGA-Kraken-17625-Samples.csv")%>%as.data.frame()
colnames(micro)[1]<-"sampleID"
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
meta3<-meta[,c(1,3,12)]
head(meta3)
colnames(meta3)[c(2,3)]<-c("tcgaID","Types")
levels(factor(meta3$Types))
meta3$Types<-gsub("^.*-","",meta3$Types)
data<-merge(meta3,micro,by="sampleID")
data[9000:9020,1:3]
write.csv(data,"../dataset/TCGA_microreads/data/TCGAmeta_microReads.csv",row.names = F)
source("../../../Git/TuMicroAnalysor/R/requirements.R")
source("../../../Git/TuMicroAnalysor/R/setTaxonClass.R")
source("../../../Git/TuMicroAnalysor/R/taxonClassSplit.R")
source("../../../Git/TuMicroAnalysor/R/reads2abundance.R")
source("../../../Git/TuMicroAnalysor/R/OTUtable_merge.R")
OTUtable<-fread("../dataset/TCGA_microreads/userProcessedData/Kraken-TCGA-OTUtable-17625-Samples.csv")%>%as.data.frame()
OTUtable[1:3,1:2]
sumreads<-data.frame(sampleID=colnames(OTUtable)[-1],sumreads=apply(OTUtable[,-1], 2, sum))
meta_sumreads<-merge(sumreads,meta1,by="sampleID")
colnames(meta_sumreads)[13]<-"tcgaID"
meta_sumreads$OS<-meta_sumreads$vital_status
meta_sumreads$OS<-ifelse(meta_sumreads$OS=="Alive",0,ifelse(meta_sumreads$OS=="Dead",1,"NA"))
meta_sumreads$OStime<-round(meta_sumreads$days_to_death/365,0)
meta_sumreads$OS<-as.numeric(meta_sumreads$OS)
meta_sumreads$OStime<-as.numeric(meta_sumreads$OStime)
write.csv(meta_sumreads,"../dataset/TCGA_microreads/userProcessedData/meta_sumreads_TCGA_10229.csv",row.names = F)
library(survminer)
library(survival)
library(survivalROC)
library(forestplot)
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
for (i in seq_along(datalist)) {
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
  pdf(file = paste0("../dataset/TCGA_microreads/userProcessedData/",names(datalist)[i],"_OS.pdf"),width = 5,height = 6,onefile = FALSE)
  print(p[[i]])
  dev.off()
  png(file = paste0("../dataset/TCGA_microreads/userProcessedData/",names(datalist)[i],"_OS.png"))
  print(p[[i]])
  dev.off()
}
p[[1]]
microScore<-fread("../dataset/TCGA_data/MRscore_TCGA_Patients9359.csv")%>%as.data.frame()
meta_sumreads$tcgaID[1]
colnames(microScore)[1]<-"tcgaID"
colnames(microScore)[1]
colnames(meta_sumreads)[4]
MIRscore_microReads<-merge(microScore,meta_sumreads[,c(1,2,4)],by="tcgaID")
write.csv(MIRscore_microReads,"../dataset/TCGA_microreads/userProcessedData/MIRscore_microReads.csv",row.names = F)
library(ggthemes)
library(ggsci)
ggplot(MIRscore_microReads,aes(MRscore,log10(sumreads)))+
  geom_point(aes(color=Types))+
  theme_few(base_size = 12)+
  scale_color_futurama()
ggplot(MIRscore_microReads,aes(MRscore,log10(sumreads)))+
  geom_point(aes(color=Types))+
  theme_bw(base_size = 12)+
  scale_color_uchicago()+
  facet_wrap(~Types,scales = "free")
data = MIRscore_microReads
x = "MRscore"
y = "sumreads"
inputFacet = "Types"
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
MIRscore_microReads$sumreads1<-log10(MIRscore_microReads$sumreads)
p<-plot_cor_with_label(data = MIRscore_microReads,
                       x = "MRscore",
                       y = "sumreads1",
                       inputFacet = "Types")
p


cor_cross_matix<-function(data,x,y,inputFacet){
  cor.coef = TRUE
  corr_eqn <- function(x,y, method='pearson', digits = 2) {
    corr_coef <- round(cor.test(x, y, method=method)$estimate, digits = digits)
    corr_pval <- tryCatch(format(cor.test(x,y, method=method)$p.value, 
                                 scientific=TRUE),
                          error=function(e) NA)
    paste(method, 'r = ', corr_coef, ',', 'p =', round(as.numeric(corr_pval),4))
  }
  if (cor.coef) {
    resCor <- data.frame(facets = unique(data[, inputFacet]))
    for(i in seq_along(resCor$facets)) {
      foo <- data[data[,inputFacet] == resCor$facets[i], ]
      resCor$text[i] <- corr_eqn(foo[,x], foo[,y])
    }}
  colnames(resCor)[1] <- inputFacet
  resCor$text<-gsub("pearson ","",resCor$text)
  resCor1<-resCor
  resCor2<-resCor
  resCor1$text<-gsub(" , p.*$","",resCor1$text)
  resCor1$text<-gsub("^.*= ","",resCor1$text)%>%as.numeric()
  resCor2$text<-gsub("^.*p","p",resCor2$text)
  resCor2$text<-gsub("^.*= ","",resCor2$text)%>%as.numeric()
  cor_res<-data.frame(CancerTypes=resCor1$Types,r_value=resCor1$text,p_value=resCor2$text)
  colnames(cor_res)[2]<-paste0(y,"_r_value")
  colnames(cor_res)[3]<-paste0(y,"_p_value")
  return(corMat=cor_res)
}


which(duplicated(meta_sumreads$investigation))
colnames(meta_sumreads)
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
genus<-fread("../dataset/TCGA_microreads/userProcessedData/Genus_abundance.csv")%>%as.data.frame()
genus[1:3,1:3]
genus<-data.frame(Genus=genus$Genus,round(genus[,-1],2))
genus<-genus%>%mutate(total_reads=apply(genus[,-1], 1,sum))
which(genus$total_reads==0)
genus<-genus[-which(genus$total_reads==0),]#1993 move 639 sum=0,remained 354 genus
grep("total_reads",colnames(genus))
genus<-genus[,c(10232,1:10231)]
genus[1:3,1:3]
both_genus<-data.frame(Genus=genus[,2],genus[,colnames(genus)%in%bothmeta$sampleID])
both_genus[1:3,1:3]
rownames(both_genus)<-both_genus$Genus
both_genus<-data.frame(t(both_genus))[-1,]
both_genus<-data.frame(sampleID=rownames(both_genus),both_genus)
both_genus<-merge(bothmeta,both_genus,by="sampleID")
write.csv(both_genus,"../dataset/TCGA_microreads/userProcessedData/MGS_RNAseq_overlap_microReads.csv")
levels(factor(both_genus$investigation))
colnames(both_genus)[41:43]
dna<-subset(both_genus,experimental_strategy=="WGS")
rna<-subset(both_genus,experimental_strategy=="RNA-Seq")
dna_meta<-data.frame(Samples=dna$filename,Group=dna$investigation)
rna_meta<-data.frame(Samples=rna$filename,Group=rna$investigation)
mciroreads<-both_genus[,c(1,43:396)]
dna_micro<-data.frame(row.names =dna$sampleID,as.matrix(dna[,c(43:396)]))
rna_micro<-data.frame(row.names =rna$sampleID,as.matrix(rna[,c(43:396)]))
library(pheatmap)
library(reshape2)
data1<-melt(rna[,c(1,12,43:396)],id.vars = c("sampleID","investigation"),variable.name = "Genus")
data2<-melt(dna[,c(1,43,12:396)],id.vars = c("sampleID","investigation"),variable.name = "Genus")
head(data1)
p1<-ggplot(data1,aes(reorder(sampleID,value),value,fill=Genus))+
  geom_col()+
  guides(fill=FALSE)+
  facet_wrap(~investigation,scales = "free")
    pheatmap(dna_micro[1:3,1:3])
p2<-ggplot(data2,aes(reorder(sampleID,value),value,fill=Genus))+
      geom_col()+
      guides(fill=FALSE)+
  facet_wrap(~investigation,scales = "free")
library(ggpubr)
ggarrange(p1,p2,ncol = 1,nrow =2)
  
data<-both_genus[,-c(2,)]

data<-data.frame(row.names = genus$Genus,genus[,-c(1:2)])
micro<-data.frame(t(data))
micro$sampleID=rownames(micro)
grep("sampleID",colnames(micro))
micro<-micro[,c(355,1:354)]
micro[1:3,1:3]
genus<-micro
meta<-fread("../dataset/TCGA_microreads/data/TCGAmeta_microReads.csv")%>%as.data.frame()
meta[1:3,1:4]
genus_meta<-merge(meta[,c(1:3)],genus,by="sampleID")
genus_meta[1:10,3:6]
write.csv(genus_meta,"../dataset/TCGA_microreads/userProcessedData/genus_meta.csv",row.names = F)
genus_meta<-fread("../dataset/TCGA_microreads/userProcessedData/genus_meta.csv")%>%as.data.frame()
sumsamples<-genus_meta[,-c(1:2)]%>%group_by(Types)%>%summarise(n())
meanreads<-data.frame(Types=genus_meta$Types,mean=apply(genus_meta[,-c(1:3)], 1, mean))
MRexpr<-read.csv("../dataset/TCGA_data/expr_MicroSignature_TMP.csv",header = T,row.names = 1)%>%as.data.frame()
MRexpr<-data.frame(t(MRexpr))
MRexpr[1:3,1:3]
MRexpr<-data.frame(sampleID=rownames(MRexpr),MRexpr)
meta<-fread("../dataset/TCGA_microreads/data/Metadata-TCGA-Kraken-17625-Samples.csv")%>%as.data.frame()
colnames(meta)[1]<-"sampleID"
meta[1:3,1:3]
meta1<-meta[grep("TCGA",meta$filename),]
meta1$filename<-gsub("^.*TCGA","TCGA",meta1$filename)
meta1$filename<-gsub("-","\\.",meta1$filename)
meta1$filename<-str_sub(meta1$filename,1,28)
meta1[1:3,8:15]
colnames(MRexpr)[1]<-"filename"
intersect(meta1$filename,MRexpr$filename)
MRsxpr_meta<-merge(meta1[,c(1,3,12)],MRexpr,by="filename")
MRsxpr_meta[1:3,1:3]
colnames(MRsxpr_meta)[3]<-"Types"
micro<-fread("../dataset/TCGA_microreads/userProcessedData/Kraken-TCGA-OTUsplitReads-17625-Samples")
genus_reads<-data.frame(row.names = micro$Genus,micro[,-c(1:9)])
genus_reads<-data.frame(t(genus_reads))
genus_reads<-data.frame(sampleID=rownames(genus_reads),genus_reads)
genus_reads[1:3,1:3]
MRsxpr_meta[1:3,1:5]
MRsxpr_genus<-merge(MRsxpr_meta,genus_reads,by="sampleID")
MRsxpr_genus$Types<-gsub("TCGA-","",MRsxpr_genus$Types)
MRsxpr_genus[1:3,1:3]
write.csv(MRsxpr_genus,"../dataset/TCGA_microreads/userProcessedData/MRexpr_genus_reads.csv")
source("../../../Git/TuMicroAnalysor/R/requirements.R")

MRexpr_genus<-fread("../dataset/TCGA_microreads/userProcessedData/MRexpr_genus_reads.csv")[,-c(1,3)]
MRexpr_genus[1:3,1:5]
data1<-split.data.frame(MRexpr_genus,f = MRexpr_genus$Types,drop = F)
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
colnames(MRexpr_genus)[3:830]
gene<-colnames(MRexpr_genus)[3:830]
micro<-colnames(MRexpr_genus)[831:2822]
library(Hmisc)
d1<-vector(mode = "list")
res<-vector(mode = "list")
cor<-vector(mode = "list")
micro_tcga_cor<-vector(mode = "list")
micro_tcga_cor1<-vector(mode = "list")
path<-vector(mode = "list")
for (i in 1:5) {
  d1[[i]]<-data1[[i]][,-c(1,2)]
  res[[i]] <- rcorr(as.matrix(d1[[i]]))
  cor[[i]]<-flattenCorrMatrix(res[[i]]$r, res[[i]]$P)
  micro_tcga_cor[[i]]<-cor[[i]] %>% filter(abs(cor)>0.7&p<0.01) %>% .[which(.$row %in% gene),] %>% .[which(.$column %in% micro),] 
  micro_tcga_cor1[[i]] <- micro_tcga_cor[[i]] %>% group_by(row) %>% mutate(GeneDegree = n()) %>% group_by(column) %>% mutate(TaxonDegree= n()) 
  path[[i]]<-paste0(names(data1)[i],"_micro_tcga_cor.csv")
  write_csv(micro_tcga_cor1[[i]],path[[i]])
}
  MRexpr<-read.csv("../dataset/TCGA_data/expr_MicroSignature_TMP.csv",header = T,row.names = 1)%>%as.data.frame()
MRexpr<-data.frame(Gene=rownames(MRexpr),MRexpr)
keypath<-fread("../dataset/TCGA_data/Nod_TLR_genes_1.csv")[,c(2,4)]
keypathwaygeneExpr<-merge(keypath,MRexpr,by="Gene")
keypathwaygeneExpr<-data.frame(row.names = keypathwaygeneExpr$Gene,keypathwaygeneExpr[,-c(1,2)])
keypathwaygeneExpr[1:3,1:3]
MRexpr<-data.frame(t(keypathwaygeneExpr))
MRexpr<-data.frame(sampleID=rownames(MRexpr),MRexpr)
MRexpr[1:3,1:3]
meta<-fread("../dataset/TCGA_microreads/data/Metadata-TCGA-Kraken-17625-Samples.csv")%>%as.data.frame()
colnames(meta)[1]<-"sampleID"
meta[1:3,1:3]
meta1<-meta[grep("TCGA",meta$filename),]
meta1$filename<-gsub("^.*TCGA","TCGA",meta1$filename)
meta1$filename<-gsub("-","\\.",meta1$filename)
meta1$filename<-str_sub(meta1$filename,1,28)
meta1[1:3,8:15]
colnames(MRexpr)[1]<-"filename"
intersect(meta1$filename,MRexpr$filename)
MRsxpr_meta<-merge(meta1[,c(1,3,12)],MRexpr,by="filename")
micro<-fread("../dataset/TCGA_microreads/userProcessedData/Kraken-TCGA-OTUsplitReads-17625-Samples")
genus_reads<-data.frame(row.names = micro$Genus,micro[,-c(1:9)])
genus_reads<-data.frame(t(genus_reads))
genus_reads<-data.frame(sampleID=rownames(genus_reads),genus_reads)
genus_reads[1:3,1:3]
MRsxpr_meta[1:3,1:5]
MRsxpr_genus<-merge(MRsxpr_meta,genus_reads,by="sampleID")
write.csv(MRsxpr_genus,"../dataset/TCGA_microreads/data/MRkeygene_genus.csv",row.names = F)
MRexpr_genus<-fread("../dataset/TCGA_microreads/data/MRkeygene_genus.csv")%>%as.data.frame()
MRexpr_genus[1:3,1:5]
MRexpr_genus<-MRexpr_genus[,-2]
colnames(MRexpr_genus)[2]<-"Types"
data1<-split.data.frame(MRexpr_genus,f = MRexpr_genus$Types,drop = F)
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
colnames(MRexpr_genus)[3:106]
gene<-colnames(MRexpr_genus)[3:106]
micro<-colnames(MRexpr_genus)[107:2009]
library(Hmisc)
d1<-vector(mode = "list")
res<-vector(mode = "list")
cor<-vector(mode = "list")
micro_tcga_cor<-vector(mode = "list")
micro_tcga_cor1<-vector(mode = "list")
path<-vector(mode = "list")
for (i in 1:5) {
  d1[[i]]<-data1[[i]][,-c(1,2)]
  res[[i]] <- rcorr(as.matrix(d1[[i]]))
  cor[[i]]<-flattenCorrMatrix(res[[i]]$r, res[[i]]$P)
  micro_tcga_cor[[i]]<-cor[[i]] %>% filter(abs(cor)>0.5&p<0.05) %>% .[which(.$row %in% gene),] %>% .[which(.$column %in% micro),] 
  micro_tcga_cor1[[i]] <- micro_tcga_cor[[i]] %>% group_by(row) %>% mutate(GeneDegree = n()) %>% group_by(column) %>% mutate(TaxonDegree= n()) 
  path[[i]]<-paste0(names(data1)[i],"_key_micro_tcga_cor.csv")
  write_csv(micro_tcga_cor1[[i]],path[[i]])
}
MRexpr_genus<-fread("../dataset/TCGA_microreads/userProcessedData/MRexpr_genus_reads.csv")[,-c(1,3)]
path<-"../dataset/TCGA_microreads/Micro_genus_cor/"
files<-list.files(path = path)
filenames<-gsub("_.*","",files)
datalist<-list()
for (i in seq_along(files)) {
  datalist[[i]]<-fread(paste0(path,files[i]))%>%as.data.frame()
  datalist[[i]]<-datalist[[i]]%>%mutate(Types=rep(filenames[1],nrow(.)))
}
MRgenecor<-bind_rows(datalist)
MRgenecor1<-subset(MRgenecor,p<=0.01&abs(cor)>=0.6)
write.csv(MRgenecor1,"../dataset/TCGA_microreads/data/MRgenecor0.6(p0.01).csv",row.names = F,quote = F)
=======
library(data.table)
library(tidyverse)
library(stringr)
micro<-read.csv("../dataset/TCGA_microreads/data/Kraken-TCGA-Raw-Data-17625-Samples.csv")%>%as.data.frame()
meta<-fread("../dataset/TCGA_microreads/data/Metadata-TCGA-Kraken-17625-Samples.csv")%>%as.data.frame()
colnames(micro)[1]<-"sampleID"
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
meta3<-meta[,c(1,3,12)]
head(meta3)
colnames(meta3)[c(2,3)]<-c("tcgaID","Types")
levels(factor(meta3$Types))
meta3$Types<-gsub("^.*-","",meta3$Types)
data<-merge(meta3,micro,by="sampleID")
data[9000:9020,1:3]
write.csv(data,"../dataset/TCGA_microreads/data/TCGAmeta_microReads.csv",row.names = F)
source("../../../Git/TuMicroAnalysor/R/requirements.R")
source("../../../Git/TuMicroAnalysor/R/setTaxonClass.R")
source("../../../Git/TuMicroAnalysor/R/taxonClassSplit.R")
source("../../../Git/TuMicroAnalysor/R/reads2abundance.R")
source("../../../Git/TuMicroAnalysor/R/OTUtable_merge.R")
OTUtable<-fread("../dataset/TCGA_microreads/userProcessedData/Kraken-TCGA-OTUtable-17625-Samples.csv")%>%as.data.frame()
OTUtable[1:3,1:2]
sumreads<-data.frame(sampleID=colnames(OTUtable)[-1],sumreads=apply(OTUtable[,-1], 2, sum))
meta_sumreads<-merge(sumreads,meta1,by="sampleID")
colnames(meta_sumreads)[13]<-"tcgaID"
meta_sumreads$OS<-meta_sumreads$vital_status
meta_sumreads$OS<-ifelse(meta_sumreads$OS=="Alive",0,ifelse(meta_sumreads$OS=="Dead",1,"NA"))
meta_sumreads$OStime<-round(meta_sumreads$days_to_death/365,0)
meta_sumreads$OS<-as.numeric(meta_sumreads$OS)
meta_sumreads$OStime<-as.numeric(meta_sumreads$OStime)
write.csv(meta_sumreads,"../dataset/TCGA_microreads/userProcessedData/meta_sumreads_TCGA_10229.csv",row.names = F)
library(survminer)
library(survival)
library(survivalROC)
library(forestplot)
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
for (i in seq_along(datalist)) {
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
  pdf(file = paste0("../dataset/TCGA_microreads/userProcessedData/",names(datalist)[i],"_OS.pdf"),width = 5,height = 6,onefile = FALSE)
  print(p[[i]])
  dev.off()
  png(file = paste0("../dataset/TCGA_microreads/userProcessedData/",names(datalist)[i],"_OS.png"))
  print(p[[i]])
  dev.off()
}
p[[1]]
microScore<-fread("../dataset/TCGA_data/MRscore_TCGA_Patients9359.csv")%>%as.data.frame()
meta_sumreads$tcgaID[1]
colnames(microScore)[1]<-"tcgaID"
colnames(microScore)[1]
colnames(meta_sumreads)[4]
MIRscore_microReads<-merge(microScore,meta_sumreads[,c(1,2,4)],by="tcgaID")
write.csv(MIRscore_microReads,"../dataset/TCGA_microreads/userProcessedData/MIRscore_microReads.csv",row.names = F)
library(ggthemes)
library(ggsci)
ggplot(MIRscore_microReads,aes(MRscore,log10(sumreads)))+
  geom_point(aes(color=Types))+
  theme_few(base_size = 12)+
  scale_color_futurama()
ggplot(MIRscore_microReads,aes(MRscore,log10(sumreads)))+
  geom_point(aes(color=Types))+
  theme_bw(base_size = 12)+
  scale_color_uchicago()+
  facet_wrap(~Types,scales = "free")
data = MIRscore_microReads
x = "MRscore"
y = "sumreads"
inputFacet = "Types"
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
MIRscore_microReads$sumreads1<-log10(MIRscore_microReads$sumreads)
p<-plot_cor_with_label(data = MIRscore_microReads,
                       x = "MRscore",
                       y = "sumreads1",
                       inputFacet = "Types")
p


cor_cross_matix<-function(data,x,y,inputFacet){
  cor.coef = TRUE
  corr_eqn <- function(x,y, method='pearson', digits = 2) {
    corr_coef <- round(cor.test(x, y, method=method)$estimate, digits = digits)
    corr_pval <- tryCatch(format(cor.test(x,y, method=method)$p.value, 
                                 scientific=TRUE),
                          error=function(e) NA)
    paste(method, 'r = ', corr_coef, ',', 'p =', round(as.numeric(corr_pval),4))
  }
  if (cor.coef) {
    resCor <- data.frame(facets = unique(data[, inputFacet]))
    for(i in seq_along(resCor$facets)) {
      foo <- data[data[,inputFacet] == resCor$facets[i], ]
      resCor$text[i] <- corr_eqn(foo[,x], foo[,y])
    }}
  colnames(resCor)[1] <- inputFacet
  resCor$text<-gsub("pearson ","",resCor$text)
  resCor1<-resCor
  resCor2<-resCor
  resCor1$text<-gsub(" , p.*$","",resCor1$text)
  resCor1$text<-gsub("^.*= ","",resCor1$text)%>%as.numeric()
  resCor2$text<-gsub("^.*p","p",resCor2$text)
  resCor2$text<-gsub("^.*= ","",resCor2$text)%>%as.numeric()
  cor_res<-data.frame(CancerTypes=resCor1$Types,r_value=resCor1$text,p_value=resCor2$text)
  colnames(cor_res)[2]<-paste0(y,"_r_value")
  colnames(cor_res)[3]<-paste0(y,"_p_value")
  return(corMat=cor_res)
}


which(duplicated(meta_sumreads$investigation))
colnames(meta_sumreads)
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
genus<-fread("../dataset/TCGA_microreads/userProcessedData/Genus_abundance.csv")%>%as.data.frame()
genus[1:3,1:3]
genus<-data.frame(Genus=genus$Genus,round(genus[,-1],2))
genus<-genus%>%mutate(total_reads=apply(genus[,-1], 1,sum))
which(genus$total_reads==0)
genus<-genus[-which(genus$total_reads==0),]#1993 move 639 sum=0,remained 354 genus
grep("total_reads",colnames(genus))
genus<-genus[,c(10232,1:10231)]
genus[1:3,1:3]
both_genus<-data.frame(Genus=genus[,2],genus[,colnames(genus)%in%bothmeta$sampleID])
both_genus[1:3,1:3]
rownames(both_genus)<-both_genus$Genus
both_genus<-data.frame(t(both_genus))[-1,]
both_genus<-data.frame(sampleID=rownames(both_genus),both_genus)
both_genus<-merge(bothmeta,both_genus,by="sampleID")
write.csv(both_genus,"../dataset/TCGA_microreads/userProcessedData/MGS_RNAseq_overlap_microReads.csv")
levels(factor(both_genus$investigation))
colnames(both_genus)[41:43]
dna<-subset(both_genus,experimental_strategy=="WGS")
rna<-subset(both_genus,experimental_strategy=="RNA-Seq")
dna_meta<-data.frame(Samples=dna$filename,Group=dna$investigation)
rna_meta<-data.frame(Samples=rna$filename,Group=rna$investigation)
mciroreads<-both_genus[,c(1,43:396)]
dna_micro<-data.frame(row.names =dna$sampleID,as.matrix(dna[,c(43:396)]))
rna_micro<-data.frame(row.names =rna$sampleID,as.matrix(rna[,c(43:396)]))
library(pheatmap)
library(reshape2)
data1<-melt(rna[,c(1,12,43:396)],id.vars = c("sampleID","investigation"),variable.name = "Genus")
data2<-melt(dna[,c(1,43,12:396)],id.vars = c("sampleID","investigation"),variable.name = "Genus")
head(data1)
p1<-ggplot(data1,aes(reorder(sampleID,value),value,fill=Genus))+
  geom_col()+
  guides(fill=FALSE)+
  facet_wrap(~investigation,scales = "free")
    pheatmap(dna_micro[1:3,1:3])
p2<-ggplot(data2,aes(reorder(sampleID,value),value,fill=Genus))+
      geom_col()+
      guides(fill=FALSE)+
  facet_wrap(~investigation,scales = "free")
library(ggpubr)
ggarrange(p1,p2,ncol = 1,nrow =2)
  
data<-both_genus[,-c(2,)]

data<-data.frame(row.names = genus$Genus,genus[,-c(1:2)])
micro<-data.frame(t(data))
micro$sampleID=rownames(micro)
grep("sampleID",colnames(micro))
micro<-micro[,c(355,1:354)]
micro[1:3,1:3]
genus<-micro
meta<-fread("../dataset/TCGA_microreads/data/TCGAmeta_microReads.csv")%>%as.data.frame()
meta[1:3,1:4]
genus_meta<-merge(meta[,c(1:3)],genus,by="sampleID")
genus_meta[1:10,3:6]
write.csv(genus_meta,"../dataset/TCGA_microreads/userProcessedData/genus_meta.csv",row.names = F)
genus_meta<-fread("../dataset/TCGA_microreads/userProcessedData/genus_meta.csv")%>%as.data.frame()
sumsamples<-genus_meta[,-c(1:2)]%>%group_by(Types)%>%summarise(n())
meanreads<-data.frame(Types=genus_meta$Types,mean=apply(genus_meta[,-c(1:3)], 1, mean))
MRexpr<-read.csv("../dataset/TCGA_data/expr_MicroSignature_TMP.csv",header = T,row.names = 1)%>%as.data.frame()
MRexpr<-data.frame(t(MRexpr))
MRexpr[1:3,1:3]
MRexpr<-data.frame(sampleID=rownames(MRexpr),MRexpr)
meta<-fread("../dataset/TCGA_microreads/data/Metadata-TCGA-Kraken-17625-Samples.csv")%>%as.data.frame()
colnames(meta)[1]<-"sampleID"
meta[1:3,1:3]
meta1<-meta[grep("TCGA",meta$filename),]
meta1$filename<-gsub("^.*TCGA","TCGA",meta1$filename)
meta1$filename<-gsub("-","\\.",meta1$filename)
meta1$filename<-str_sub(meta1$filename,1,28)
meta1[1:3,8:15]
colnames(MRexpr)[1]<-"filename"
intersect(meta1$filename,MRexpr$filename)
MRsxpr_meta<-merge(meta1[,c(1,3,12)],MRexpr,by="filename")
MRsxpr_meta[1:3,1:3]
colnames(MRsxpr_meta)[3]<-"Types"
micro<-fread("../dataset/TCGA_microreads/userProcessedData/Kraken-TCGA-OTUsplitReads-17625-Samples")
genus_reads<-data.frame(row.names = micro$Genus,micro[,-c(1:9)])
genus_reads<-data.frame(t(genus_reads))
genus_reads<-data.frame(sampleID=rownames(genus_reads),genus_reads)
genus_reads[1:3,1:3]
MRsxpr_meta[1:3,1:5]
MRsxpr_genus<-merge(MRsxpr_meta,genus_reads,by="sampleID")
MRsxpr_genus$Types<-gsub("TCGA-","",MRsxpr_genus$Types)
MRsxpr_genus[1:3,1:3]
write.csv(MRsxpr_genus,"../dataset/TCGA_microreads/userProcessedData/MRexpr_genus_reads.csv")
source("../../../Git/TuMicroAnalysor/R/requirements.R")

MRexpr_genus<-fread("../dataset/TCGA_microreads/userProcessedData/MRexpr_genus_reads.csv")[,-c(1,3)]
MRexpr_genus[1:3,1:5]
data1<-split.data.frame(MRexpr_genus,f = MRexpr_genus$Types,drop = F)
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
colnames(MRexpr_genus)[3:830]
gene<-colnames(MRexpr_genus)[3:830]
micro<-colnames(MRexpr_genus)[831:2822]
library(Hmisc)
d1<-vector(mode = "list")
res<-vector(mode = "list")
cor<-vector(mode = "list")
micro_tcga_cor<-vector(mode = "list")
micro_tcga_cor1<-vector(mode = "list")
path<-vector(mode = "list")
for (i in 1:5) {
  d1[[i]]<-data1[[i]][,-c(1,2)]
  res[[i]] <- rcorr(as.matrix(d1[[i]]))
  cor[[i]]<-flattenCorrMatrix(res[[i]]$r, res[[i]]$P)
  micro_tcga_cor[[i]]<-cor[[i]] %>% filter(abs(cor)>0.7&p<0.01) %>% .[which(.$row %in% gene),] %>% .[which(.$column %in% micro),] 
  micro_tcga_cor1[[i]] <- micro_tcga_cor[[i]] %>% group_by(row) %>% mutate(GeneDegree = n()) %>% group_by(column) %>% mutate(TaxonDegree= n()) 
  path[[i]]<-paste0(names(data1)[i],"_micro_tcga_cor.csv")
  write_csv(micro_tcga_cor1[[i]],path[[i]])
}
  MRexpr<-read.csv("../dataset/TCGA_data/expr_MicroSignature_TMP.csv",header = T,row.names = 1)%>%as.data.frame()
MRexpr<-data.frame(Gene=rownames(MRexpr),MRexpr)
keypath<-fread("../dataset/TCGA_data/Nod_TLR_genes_1.csv")[,c(2,4)]
keypathwaygeneExpr<-merge(keypath,MRexpr,by="Gene")
keypathwaygeneExpr<-data.frame(row.names = keypathwaygeneExpr$Gene,keypathwaygeneExpr[,-c(1,2)])
keypathwaygeneExpr[1:3,1:3]
MRexpr<-data.frame(t(keypathwaygeneExpr))
MRexpr<-data.frame(sampleID=rownames(MRexpr),MRexpr)
MRexpr[1:3,1:3]
meta<-fread("../dataset/TCGA_microreads/data/Metadata-TCGA-Kraken-17625-Samples.csv")%>%as.data.frame()
colnames(meta)[1]<-"sampleID"
meta[1:3,1:3]
meta1<-meta[grep("TCGA",meta$filename),]
meta1$filename<-gsub("^.*TCGA","TCGA",meta1$filename)
meta1$filename<-gsub("-","\\.",meta1$filename)
meta1$filename<-str_sub(meta1$filename,1,28)
meta1[1:3,8:15]
colnames(MRexpr)[1]<-"filename"
intersect(meta1$filename,MRexpr$filename)
MRsxpr_meta<-merge(meta1[,c(1,3,12)],MRexpr,by="filename")
micro<-fread("../dataset/TCGA_microreads/userProcessedData/Kraken-TCGA-OTUsplitReads-17625-Samples")
genus_reads<-data.frame(row.names = micro$Genus,micro[,-c(1:9)])
genus_reads<-data.frame(t(genus_reads))
genus_reads<-data.frame(sampleID=rownames(genus_reads),genus_reads)
genus_reads[1:3,1:3]
MRsxpr_meta[1:3,1:5]
MRsxpr_genus<-merge(MRsxpr_meta,genus_reads,by="sampleID")
write.csv(MRsxpr_genus,"../dataset/TCGA_microreads/data/MRkeygene_genus.csv",row.names = F)
MRexpr_genus<-fread("../dataset/TCGA_microreads/data/MRkeygene_genus.csv")%>%as.data.frame()
MRexpr_genus[1:3,1:5]
MRexpr_genus<-MRexpr_genus[,-2]
colnames(MRexpr_genus)[2]<-"Types"
data1<-split.data.frame(MRexpr_genus,f = MRexpr_genus$Types,drop = F)
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
colnames(MRexpr_genus)[3:106]
gene<-colnames(MRexpr_genus)[3:106]
micro<-colnames(MRexpr_genus)[107:2009]
library(Hmisc)
d1<-vector(mode = "list")
res<-vector(mode = "list")
cor<-vector(mode = "list")
micro_tcga_cor<-vector(mode = "list")
micro_tcga_cor1<-vector(mode = "list")
path<-vector(mode = "list")
for (i in 1:5) {
  d1[[i]]<-data1[[i]][,-c(1,2)]
  res[[i]] <- rcorr(as.matrix(d1[[i]]))
  cor[[i]]<-flattenCorrMatrix(res[[i]]$r, res[[i]]$P)
  micro_tcga_cor[[i]]<-cor[[i]] %>% filter(abs(cor)>0.5&p<0.05) %>% .[which(.$row %in% gene),] %>% .[which(.$column %in% micro),] 
  micro_tcga_cor1[[i]] <- micro_tcga_cor[[i]] %>% group_by(row) %>% mutate(GeneDegree = n()) %>% group_by(column) %>% mutate(TaxonDegree= n()) 
  path[[i]]<-paste0(names(data1)[i],"_key_micro_tcga_cor.csv")
  write_csv(micro_tcga_cor1[[i]],path[[i]])
}
MRexpr_genus<-fread("../dataset/TCGA_microreads/userProcessedData/MRexpr_genus_reads.csv")[,-c(1,3)]
path<-"../dataset/TCGA_microreads/Micro_genus_cor/"
files<-list.files(path = path)
filenames<-gsub("_.*","",files)
datalist<-list()
for (i in seq_along(files)) {
  datalist[[i]]<-fread(paste0(path,files[i]))%>%as.data.frame()
  datalist[[i]]<-datalist[[i]]%>%mutate(Types=rep(filenames[1],nrow(.)))
}
MRgenecor<-bind_rows(datalist)
MRgenecor1<-subset(MRgenecor,p<=0.01&abs(cor)>=0.6)
write.csv(MRgenecor1,"../dataset/TCGA_microreads/data/MRgenecor0.6(p0.01).csv",row.names = F,quote = F)
>>>>>>> 8bb5eb81eae04d21b86a863239c44287b94d2c46
    