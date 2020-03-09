library(tidyverse)
library(data.table)
library(GEOquery)
options(stringsAsFactors = F)
#mc1="GSE13015.Rdata"
mc3="GSE40012.Rdata"
mc5="GSE65682.Rdata"
mc6="GSE21802.Rdata"
#mc7="GSE27131.Rdata"#few data,1456 genes
mc8="GSE28750.Rdata"
mc9="GSE42834.Rdata"
mc10="GSE57065.Rdata"
mc11="GSE68310.Rdata"
mc12="GSE69528.Rdata"
mc13="GSE82050.Rdata"
mc14="GSE111368.Rdata"
datalist<-list(mc3,mc5,mc6,mc8,mc10,mc11,mc12,mc13,mc14)
GEOid<-gsub(".Rdata","",datalist)
downloaded<-list.files("../dataset/dataset_alidation/validate_datasets/")
downloadedID<-levels(factor(gsub("_.*","",downloaded)))
which(GEOid%in%downloadedID)
expr_list<-paste0("../dataset/dataset_alidation/validate_datasets/",GEOid,"_expr.csv")
meta_list<-paste0("../dataset/dataset_alidation/validate_datasets/",GEOid,"_meta.csv")
head(expr_list)
head(meta_list)

check_genename<-function(expr){
  print(expr[1:6,1])
}
for (i in seq(length(expr_list))) {
  print(i)
  check_genename(fread(expr_list[i]))
}

##GEOid[c(1,3,4,5,6,9,10)] have null value in the first row
deal_id<-expr_list[c(1,3,4,5,6,9,10)]
deal_expr<-list()
for (i in seq(length(deal_id))) {
  print(i)
  deal_expr[[i]]<-fread(deal_id[i])
  deal_expr[[i]]<-deal_expr[[i]][-1,]
  write.csv(deal_expr[[i]],deal_id[i],row.names = T)
}
for (i in seq(length(expr_list))) {
  print(i)
  check_genename(fread(expr_list[i]))
}
expr_list<-paste0("../dataset/dataset_alidation/validate_datasets/",GEOid,"_expr.csv")
exprall<-list()
metall<-list()
for (i in seq(length(expr_list))) {
  print(i)
  exprall[[i]]<-fread(expr_list[i])
  metall[[i]]<-fread(meta_list[i])
  print(exprall[[i]][1:10,1:3])
  names(exprall)[i]<-GEOid[i]
  names(metall)[i]<-GEOid[i]
}
#exprall[[c(1,3)]] haven't been normalized
logNormal<-function(x){
  log(x+1)
}
exprall[[1]]<-data.frame(Symbol=exprall[[1]][,1],apply(exprall[[1]][,-1], 2, logNormal))
exprall[[3]]<-data.frame(Symbol=exprall[[3]][,1],apply(exprall[[3]][,-1], 2, logNormal))
sapply(exprall, dim)
#       GSE40012 GSE65682 GSE21802 GSE28750 GSE57065 GSE68310 GSE69528 GSE82050
#[1,]    25137     1909    24614    23347    23348      579     1712     2317
#[2,]      191      803       41       42      108      881      139       40
names(exprall)[1]
names(metall)[1]

largeExpr<-inner_join(inner_join(inner_join(exprall[[1]],exprall[[3]]),exprall[[4]]),exprall[[5]])
lessExpr<-inner_join(inner_join(exprall[[2]],exprall[[8]]),exprall[[9]])
View(metall[[1]])
group1<-data.frame(sampleID=metall[[1]]$geo_accession,
                   Day=metall[[1]]$`day:ch1`,
                   Age=metall[[1]]$`ID:ch1`,
                   Platform_id=metall[[1]]$platform_id,
                   Gender=metall[[1]]$`gender:ch1`,
                   Group=metall[[1]]$`sample type:ch1`,
                   Tissue=metall[[1]]$`tissue:ch1`)
group1$Tissue<-rep("Whole blood",nrow(group1))
group1$Type<-group1$Group
levels(factor(group1$Type))
group1$Type<-gsub("influenza A pneumonia","Viral_infection",group1$Type)
group1$Type<-gsub("healthy control","Healthy",group1$Type)
group1$Type<-gsub("bacterial pneumonia","Bacterial_infection",group1$Type)
group1$Type<-gsub("mixed bacterial and influenza A pneumonia","Mixed_infection",group1$Type)
group1$Type<-gsub("systemic inflammatory response","Uninfected_inflammation",group1$Type)
group1$Type<-gsub("mixed bacterial and Viral_infection","Mixed_infection",group1$Type)
group1$Healthy0.Infection<-group1$Type
group1$Healthy0.Infection<-ifelse(group1$Type=="Healthy",0,1)
head(group1)

View(metall[[2]])
names(metall[2])
group2<-data.frame(sampleID=metall[[2]]$geo_accession,
                   Day=metall[[2]]$`time_to_event_28days:ch1`,
                   Platform_id=metall[[2]]$platform_id,
                   Age=metall[[2]]$`age:ch1`,
                   Gender=metall[[2]]$`gender:ch1`,
                   Group=metall[[2]]$characteristics_ch1.2,
                   Group1=metall[[2]]$`icu_acquired_infection:ch1`,
                   Group2=metall[[2]]$`endotype_class:ch1`,
                   Tissue=metall[[2]]$source_name_ch1)
group2$Type<-group2$Group
levels(factor(group2$Group))
group2$Type<-ifelse(group2$Type=="pneumonia diagnoses: NA","Healthy","Infection")
group2$Type<-ifelse(group2$Group=="pneumonia diagnoses: no-cap","Uninfected_inflammation",group2$Type)
group2$Type<-ifelse(group2$Group=="pneumonia diagnoses: cap"|group2$Group=="pneumonia diagnoses: hap","Bacterial_infection",group2$Type)
group2$Healthy0.Infection<-group2$Type
group2$Healthy0.Infection<-ifelse(group2$Type=="Healthy",0,1)
head(group2)
View(metall[[3]])
group3<-data.frame(sampleID=metall[[3]]$geo_accession,
                   Group=metall[[3]]$`virus strain:ch1`,
                   Platform_id=metall[[3]]$platform_id,
                   Group1=metall[[3]]$`disease phase:ch1`,
                   Tissue=metall[[3]]$`tissue:ch1`)
group3$Type<-group3$Group
group3$Type<-ifelse(group3$Type=="Pandemic influenza (A/H1N1 new subtype)","Viral_infection","Healthy")
group3$Healthy0.Infection<-group3$Type
group3$Healthy0.Infection<-ifelse(group3$Type=="Healthy",0,1)
head(group3)
names(metall[4])
View(metall[[4]])
group4<-data.frame(sampleID=metall[[4]]$geo_accession,
                   Group=metall[[4]]$`health status:ch1`,
                   Platform_id=metall[[4]]$platform_id,
                   Tissue=metall[[4]]$`tissue:ch1`)
group4$Type<-group4$Group
levels(factor(group4$Type))
group4$Type<-gsub("HEALTHY","Healthy",group4$Type)
group4$Type<-gsub("SEPSIS","Sepsis",group4$Type)
group4$Type<-gsub("POST_SURGICAL","Sepsis_surgical",group4$Type)
group4$Healthy0.Infection<-group4$Type
group4$Healthy0.Infection<-ifelse(group4$Type=="Healthy",0,1)
head(group4)
View(metall[[5]])
#Early and dynamic changes in gene expression in septic shock 败血性休克 patients
#28 ICU patients were enrolled at the onset of septic shock.
#Blood samples were collected within 30 minutes, 24 and 48 hours after shock 
# and compared to 25 healthy volunteers
group5<-data.frame(sampleID=metall[[5]]$geo_accession,
                   Group=metall[[5]]$`sapsii:ch1`,
                   Time=metall[[5]]$`collection time:ch1`,
                   Platform_id=metall[[5]]$platform_id,
                   Age=metall[[5]]$`age:ch1`,
                   Gender=metall[[5]]$`gender:ch1`,
                   Tissue=metall[[5]]$extract_protocol_ch1)
group5$Tissue<-ifelse(group5$Tissue=="Whole blood Paxgene tubes extraction","Whole blood",group5$Tissue)
group5$Type<-group5$Group
group5$Type<-ifelse(group5$Type=="NA","Healthy","Sepsis")
group5$Healthy0.Infection<-group5$Type
group5$Healthy0.Infection<-ifelse(group5$Type=="Healthy",0,1)
head(group5)
View(metall[[6]])
#	Healthy adults were invited to enroll to be followed for acute respiratory illness (ARI) 
#Seventy-three (55%) had an influenza virus infection, 
#64 influenza A and 9 influenza B. 
#The remaining subjects had a rhinovirus infection (N = 32), 
#other viral infections (N = 4), or no viral agent identified (N = 24).
#The results, which were replicated between two seasons, 
#showed a dramatic upregulation of interferon pathway and innate immunity genes. 
#This persisted for 2-4 days. 
#The data show a recovery phase at days 4 and 6 with differentially expressed transcripts implicated in cell proliferation and repair. 
#By day 21 the gene expression pattern was indistinguishable from baseline (enrollment).
group6<-data.frame(sampleID=metall[[6]]$geo_accession,
                   Group=metall[[6]]$`infection:ch1`,
                   Platform_id=metall[[6]]$platform_id,
                   PAtientID=metall[[6]]$`subject id:ch1`,
                   Time=metall[[6]]$`time point:ch1`,
                   Gender=metall[[6]]$`gender:ch1`,
                   Tissue=metall[[6]]$source_name_ch1)
group6$Type<-group6$Group
levels(factor(group6$Type))
group6$Type<-ifelse(group6$Type=="our tests did not detect one of the viruses sought","Healthy","Viral_infection")
group6$Healthy0.Infection<-group6$Type
group6$Healthy0.Infection<-ifelse(group6$Type=="Healthy",0,1)
head(group6)
View(metall[[7]])
names(metall[7])
#29 patients with septicemic melioidosis 
#(Melioidosis, caused by Gram negative bacteria Burkholderia pseudomallei,
#is a major type of community-acquired septicemia in Southeast Asia and Northern Australia with high mortality and morbidity rate.
#and 54 patients with sepsis due to other pathogens
group7<-data.frame(sampleID=metall[[7]]$geo_accession,
                   Group=metall[[7]]$`pathogens:ch1`,
                   Platform_id=metall[[7]]$platform_id,
                   Group1=metall[[7]]$`study group:ch1`,
                   Tissue=metall[[7]]$`tissue:ch1`)
group7$Type<-group7$Group
levels(factor(group7$Type))
group7$Type<-ifelse(group7$Type=="Control","Healthy","Bacterial_infection")
group7$Healthy0.Infection<-group7$Type
group7$Healthy0.Infection<-ifelse(group7$Type=="Healthy",0,1)
head(group7)

View(metall[[8]])
names(metall[8])
group8<-data.frame(sampleID=metall[[8]]$geo_accession,
                   Group=metall[[8]]$`infection:ch1`,
                   Platform_id=metall[[8]]$platform_id,
                   Age=metall[[8]]$`age:ch1`,
                   Tissue=metall[[8]]$`tissue:ch1`)
group8$Type<-group8$Group
levels(factor(group8$Type))
group8$Type<-ifelse(group8$Type=="healthy control","Healthy","Viral_infection")
group8$Healthy0.Infection<-group8$Type
group8$Healthy0.Infection<-ifelse(group8$Type=="Healthy",0,1)
head(group8)
View(metall[[9]])
names(metall[9])
group9<-data.frame(sampleID=metall[[9]]$geo_accession,
                   Group=metall[[9]]$`flu_type:ch1`,
                   Platform_id=metall[[9]]$platform_id,
                   Group1=metall[[9]]$`bacterial_status:ch1`,
                   Age=metall[[9]]$`age:ch1`,
                   Gender=metall[[9]]$`Sex:ch1`,
                   Time=metall[[9]]$`day of illness:ch1`,
                   Tissue=metall[[9]]$source_name_ch1)
group9$Group1<-paste0("BacteriaStatus_",group9$Group1)
group9$Type<-group9$Group
levels(factor(group9$Type))
levels(factor(group9$Group1))
group9$Type<-ifelse(group9$Type=="HC","Healthy","Viral_infection")
group9$Type<-ifelse(group9$Group1=="BacteriaStatus_Yes","Mixed_infection",group9$Type)
group9$Healthy0.Infection<-group9$Type
group9$Healthy0.Infection<-ifelse(group9$Type=="Healthy",0,1)
head(group9)
rm(group)
groupList<-list(group1,group2,group3,group4,group5,group6,group7,group8,group9)
names(groupList)<-GEOid
groupNames<-paste0("../dataset/dataset_alidation/validate_datasets/",GEOid,"_group.csv")
groupNames[1]
for (i in seq(length(groupList))) {
  groupList[[i]]$Dataset<-rep(GEOid[i],nrow(groupList[[i]]))
}
for (i in seq(length(groupNames))) {
  write.csv(groupList[[i]],groupNames[i],row.names = F)
}
sapply(groupList, dim)
#  GSE40012 GSE65682 GSE21802 GSE28750 GSE57065 GSE68310 GSE69528
#[1,] 190      802       40       41      107      880      138
#[2,] 6        8        4        3        6        6        4
#GSE82050 GSE111368
#[1,]39       359
#[2,]4         7
sapply(groupList, colnames)
mergeGroup<-bind_rows(groupList[[1]][,c("sampleID","Type","Dataset","Platform_id","Healthy0.Infection")],
                      groupList[[2]][,c("sampleID","Type","Dataset","Platform_id","Healthy0.Infection")],
                      groupList[[3]][,c("sampleID","Type","Dataset","Platform_id","Healthy0.Infection")],
                      groupList[[4]][,c("sampleID","Type","Dataset","Platform_id","Healthy0.Infection")],
                      groupList[[5]][,c("sampleID","Type","Dataset","Platform_id","Healthy0.Infection")],
                      groupList[[6]][,c("sampleID","Type","Dataset","Platform_id","Healthy0.Infection")],
                      groupList[[7]][,c("sampleID","Type","Dataset","Platform_id","Healthy0.Infection")],
                      groupList[[8]][,c("sampleID","Type","Dataset","Platform_id","Healthy0.Infection")],
                      groupList[[9]][,c("sampleID","Type","Dataset","Platform_id","Healthy0.Infection")])
head(mergeGroup)
groupStat<-mergeGroup%>%group_by(Type)%>%summarise(SampleLength=n())
groupStat1<-mergeGroup%>%group_by(Type,Dataset)%>%summarise(SampleLength=n())
groupStat
groupStat1
write.csv(mergeGroup,"../dataset/dataset_alidation/validate_datasets/mergeGroup.csv",row.names = F)
write.csv(mergeGroup,"../dataset/dataset_alidation/validate_datasets/MicroInfectiondatasets_stat.csv",row.names = F)
library(ggsci)
library(ggthemes)
ggplot(groupStat,aes(reorder(Type,SampleLength),SampleLength,fill=SampleLength))+
  geom_col()+
  scale_fill_gsea()+
  theme_stata()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank())+
  geom_text(aes(x=Type,y=SampleLength+30,label=SampleLength))
ggplot(groupStat1,aes(reorder(Type,SampleLength),SampleLength,fill=SampleLength))+
  geom_col()+
  scale_fill_gsea()+
  theme_stata()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1,vjust = 1),
        axis.title.x = element_blank())+
  geom_text(aes(x=Type,y=SampleLength+30,label=SampleLength))+
  facet_grid(~Dataset,scales = "free")

##Rdata make
names(exprall)
names(groupList)
exprall[[2]][1:5,1:5]
groupList[[2]][1:5,1:5]
nameRow<-function(x){
  x<-data.frame(row.names = x$Symbol,x[,-1])
  return(x)
}
exprall<-lapply(exprall,nameRow)
exprall1[[6]][1:5,1:5]
for (i in seq(length(groupList))) {
  rownames(groupList[[i]])<-groupList[[i]]$sampleID
}
groupList[[1]][1:5,1:5]
filter_GSElist_name<-paste0("filter_",GEOid,".Rdata")
filesPath<-paste0("../dataset/dataset_alidation/filteredGSE/",filter_GSElist_name)
filter_GSElist_name
filesPath
GSEs1<-list(pheno=groupList[[1]],genes=exprall[[1]])
GSEs2<-list(pheno=groupList[[2]],genes=exprall[[2]])
GSEs3<-list(pheno=groupList[[3]],genes=exprall[[3]])
GSEs4<-list(pheno=groupList[[4]],genes=exprall[[4]])
GSEs5<-list(pheno=groupList[[5]],genes=exprall[[5]])
GSEs6<-list(pheno=groupList[[6]],genes=exprall[[6]])
GSEs7<-list(pheno=groupList[[7]],genes=exprall[[7]])
GSEs8<-list(pheno=groupList[[8]],genes=exprall[[8]])
GSEs9<-list(pheno=groupList[[9]],genes=exprall[[9]])
GSEs<-list(GSEs1,GSEs2,GSEs3,GSEs4,GSEs5,GSEs6,GSEs7,GSEs8,GSEs9)
names(GSEs)<-GEOid

save(GSEs1,file=filesPath[1])
save(GSEs2,file=filesPath[2])
save(GSEs3,file=filesPath[3])
save(GSEs4,file=filesPath[4])
save(GSEs5,file=filesPath[5])
save(GSEs6,file=filesPath[6])
save(GSEs7,file=filesPath[7])
save(GSEs8,file=filesPath[8])
save(GSEs9,file=filesPath[9])
GSEs[[1]]$pheno[1:5,1:5]
GSEs[[2]]$genes[1:5,1:5]
datasetsInfo<-list(GEOid=GEOid,filesPath=filesPath)
save(datasetsInfo,file = "../dataset/dataset_alidation/filteredGSE/datasetsInfo.Rdata")
load(file = "../dataset/dataset_alidation/filteredGSE/datasetsInfo.Rdata")
filesPath<-datasetsInfo$filesPath
datalist<-list.files("../dataset/dataset_alidation/filteredGSE/")
filesPath<-paste0("../dataset/dataset_alidation/filteredGSE/",datalist)
for (i in seq(length(filesPath))) {
  load(file = filesPath[i])
}
GSEs<-list(GSEs1,GSEs2,GSEs3,GSEs4,GSEs5,GSEs6,GSEs7,GSEs8,GSEs9)
names(GSEs)<-GEOid

##co_normalize
library(COCONUT)
sapply(GSEs, function(x){nrow(x$genes)})
#GSE40012  GSE65682  GSE21802  GSE28750  GSE57065  GSE68310  GSE69528 
#25137      1909     24614     23347     23348       579      1712 
#GSE82050 GSE111368 
#2317      1696 
largeGSEs<-GSEs[c(1,3,4)]
GSEs_COCONUT <- COCONUT(GSEs=largeGSEs,
                        control.0.col="Healthy0.Infection",
                        byPlatform=FALSE)
GSEs.COCO.combined <- combineCOCOoutput(GSEs_COCONUT)
str(GSEs.COCO.combined)
GSEs.COCO.combined$genes[1:5,1:5]
GSEs.COCO.combined$pheno[1:5,1:5]
write.csv(GSEs.COCO.combined$genes,"../dataset/dataset_alidation/validation_results/coco_3LargeExpr.csv",row.names = T)
write.csv(GSEs.COCO.combined$pheno,"../dataset/dataset_alidation/validation_results/coco_3LargeMeta.csv",row.names = T)
##plot stat of the 3 large datasets
meta<-GSEs.COCO.combined$pheno
expr<-GSEs.COCO.combined$genes
levels(factor(meta$Type))
stat_meta<-meta%>%group_by(Type)%>%summarise(SampleLength=n())
stat_meta1<-meta%>%group_by(Type,Dataset)%>%summarise(SampleLength=n())
ggplot(stat_meta,aes(reorder(Type,SampleLength),SampleLength,fill=SampleLength))+
  geom_col()+
  scale_fill_gsea()+
  theme_stata()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank())+
  geom_text(aes(x=Type,y=SampleLength+10,label=SampleLength))
ggplot(stat_meta1,aes(reorder(Type,SampleLength),SampleLength,fill=SampleLength))+
  geom_col()+
  scale_fill_gsea()+
  theme_stata()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1,vjust = 1),
        axis.title.x = element_blank())+
  geom_text(aes(x=Type,y=SampleLength+5,label=SampleLength))+
  facet_wrap(~Dataset,scales = "free")
##bar plot of the nomalized data
expr<-as.matrix(expr)
gsms <- str_c(meta$Healthy0.Infection,sep = "",collapse = "")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")
fl <- as.factor(sml)
labels <- c("Healthy","Event")
palette(c("#f4dfdf","#dfeaf4", "#AABBCC"))
pdf("../dataset/dataset_alidation/validation_results/coco-bar.pdf",width = 25,height = 6)
boxplot(expr, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")
dev.off()
##Check batch effects
library(tidyverse)
library(EnhancedVolcano)
library(factoextra)
expr<-read.csv("../dataset/dataset_alidation/validation_results/coco_3LargeExpr.csv",row.names = 1,header = T)
meta<-read.csv("../dataset/dataset_alidation/validation_results/coco_3LargeMeta.csv",row.names = 1,header = T)
expr[1:5,1:5]
meta[1:5,1:5]
expr<-as.matrix(expr)
levels(factor(meta$Type))
res.pca <- prcomp(t(expr), scale = TRUE)
PCA<-fviz_pca_ind(res.pca,
                  label = "none",
                  habillage = meta$Dataset,
                  addEllipses = TRUE,
                  ggtheme = theme_minimal(),
                  ellipse.type = "confidence")
print(PCA)
PCA<-fviz_pca_ind(res.pca,
                  label = "none",
                  habillage = meta$Type,
                  addEllipses = TRUE,
                  ggtheme = theme_minimal(),
                  ellipse.type = "confidence")
print(PCA)
##need to move batch effetcs from different datasets by sva
library(sva)
meta%>%group_by(Type,Dataset)%>%summarise(n())
meta$bachType<-sample(c(1,2,3),nrow(meta),replace = T)
meta$bachType<-ifelse(meta$Dataset=="GSE28750",1,ifelse(meta$Dataset=="GSE40012",2,3))
head(meta)
cbdata<-as.matrix(expr)
dist_mat <- dist(t(cbdata))
clustering <- hclust(dist_mat, method = "complete")
plot(clustering, labels = meta$bachType)
plot(clustering, labels =meta$Type)
mod = model.matrix(~as.factor(Dataset), data=meta)
n.sv = num.sv(cbdata,mod,method="leek")
combat_edata = ComBat(dat=cbdata, batch=meta$bachType)
par(cex = 0.7)
n.sample=ncol(combat_edata)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
boxplot(combat_edata,col=cols,main="expression value",las=2)
boxplot(combat_edata, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")
##PCA again
res.pca <- prcomp(t(combat_edata), scale = TRUE)
PCA<-fviz_pca_ind(res.pca,
                  label = "none",
                  habillage = meta$Dataset,
                  addEllipses = TRUE,
                  ggtheme = theme_minimal(),
                  ellipse.type = "confidence")
print(PCA)
PCA<-fviz_pca_ind(res.pca,
                  label = "none",
                  habillage = meta$Type,
                  addEllipses = TRUE,
                  ggtheme = theme_minimal(),
                  ellipse.type = "confidence")
print(PCA)
combat_edata[1:5,1:5]
write.csv(combat_edata,"../dataset/dataset_alidation/validation_results/coco_sva_expr.csv",row.names = T)
##DEGs by limma each-event to Healthy
options(stringsAsFactors = F)
library(limma)
library(tidyverse)
library(data.table)
library(EnhancedVolcano)
combat_edata<-read.csv("../dataset/dataset_alidation/validation_results/coco_sva_expr.csv",row.names = 1,header = T)
meta<-read.csv("../dataset/dataset_alidation/validation_results/coco_3LargeMeta.csv",row.names = 1,header = T)
levels(factor(meta$Type))
#### Viral_infection-Healthy
Type<-levels(factor(meta$Type))
Type
Type<-Type[-2]
expr2DEG<-function(event){
meta_sub<-meta[which(meta$Type%in%c("Healthy",Type[i])),]
combat_edata_sub<-combat_edata[,which(colnames(combat_edata)%in%meta_sub$sampleID)]%>%as.matrix()
design=model.matrix(~factor(meta_sub$Type)+0)
colnames(design)=levels(factor(meta_sub$Type))
mycompare<-str_c(colnames(design),collapse = "-")
contrast.matrix<-makeContrasts(mycompare,
                               levels = design)
fit=lmFit(combat_edata_sub,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2) 
DEG<-topTable(fit2, adjust="fdr",coef=1, n=Inf) %>% na.omit()  ## coef比较分组 n基因数
message(paste0("finish DEG for",event," vs Healthy analysis !"))
# Volcano plot
VolcanoPlot<-EnhancedVolcano(DEG,
                             lab = rownames(DEG),
                             x = "logFC",
                             y = "adj.P.Val",
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

pdf(paste0("../dataset/dataset_alidation/validation_results/",event,"_DEGres.pdf"))
print(VolcanoPlot)
dev.off()
png(paste0("../dataset/dataset_alidation/validation_results/",event,"_DEGres.png"))
print(VolcanoPlot)
dev.off()
write.csv(DEG,paste0("../dataset/dataset_alidation/validation_results/",event,"_DEGres.csv"),row.names = T)
return(DEG)
}

DEGres<-list()
DEGresSig<-list()
for (i in seq(length(Type))) {
  DEGres<-expr2DEG(event = Type[i])
  DEGresSig[[i]]<-subset(DEGres[[i]],adj.P.Val<=0.05)
  names(DEGres)[i]<-Type[i]
  names(DEGresSig)[i]<-Type[i]
}
DEGres<-list()
DEGresSig<-list()
for (i in seq(length(Type))) {
  DEGres[[i]]<-read.csv(paste0("../dataset/dataset_alidation/validation_results/",Type[i],"_DEGres.csv"),row.names = 1,header = T)
  DEGresSig[[i]]<-subset(DEGres[[i]],adj.P.Val<=0.05)
  DEGresSig[[i]]<-data.frame(Gene=rownames(DEGresSig[[i]]),DEGresSig[[i]])
  names(DEGres)[i]<-Type[i]
  names(DEGresSig)[i]<-Type[i]
  }
names(DEGresSig)
DEGresSig[[1]][1:5,1:5]
sapply(DEGresSig, nrow)
#Bacterial_infection         Mixed_infection 
#5913                    4867 
#Sepsis         Sepsis_surgical 
#5992                    5649 
#Uninfected_inflammation         Viral_infection 
#6037                    6380

meta_sub<-list()
myexpr<-list()
for (i in seq(length(Type))) {
  meta_sub[[i]]<-meta[which(meta$Type%in%c("Healthy",Type[i])),]
  myexpr[[i]]<-combat_edata[,which(colnames(combat_edata)%in%meta_sub[[i]]$sampleID)]%>%as.data.frame()
  myexpr[[i]]<-data.frame(Gene=rownames(myexpr[[i]]),myexpr[[i]])
  names(myexpr)[i]<-Type[i]
}
myexpr[[1]][1:5,1:5]
DEGresSig[[3]][1:5,1:5]
expr<-expr$Bacterial_infection
DEexpr<-DEGresSig$Bacterial_infection

expr2MRscore<-function(expr,DEexpr){
  signature<-fread("../dataset/TCGA_data/annotationRow1.csv")
  MRgene_DEGs<-subset(DEexpr,Gene%in%signature$Gene)
  MRgene_DEGs$Regulation<-sample(c("Up","Down"),nrow(MRgene_DEGs),replace = T)
  MRgene_DEGs$Regulation<-ifelse(MRgene_DEGs$logFC>=0,"Up","Down")
  DE_Mgene<-MRgene_DEGs[,c("Gene","Regulation")]
  Total_MRscore<-sum(subset(MRgene_DEGs,logFC>=0)$logFC)+sum(subset(MRgene_DEGs,logFC<=0)$logFC)
  message(paste0("Total Score = ",Total_MRscore))
  upgene<-subset(DE_Mgene,Regulation=="Up")
  downgene<-subset(DE_Mgene,Regulation=="Down")
  Upexpr<-expr[which(expr$Gene%in%upgene$Gene),]
  Downexpr<-expr[which(expr$Gene%in%downgene$Gene),]
  Upexpr_average<-data.frame(sampleID = colnames(Upexpr)[-1],Upscore=apply(Upexpr[,-1], 2, mean))
  Downexpr_average<-data.frame(sampleID = colnames(Downexpr)[-1],Downscore=apply(Downexpr[,-1], 2, mean))
  DE_average<-bind_cols(Upexpr_average,Downexpr_average)[,-3]
  MRscore_value<-DE_average%>%mutate(MRscore=Upscore-Downscore)
  return(list(MRgene_DEGs=MRgene_DEGs,MRscore_value=MRscore_value,Total_MRscore=Total_MRscore))
}
names(myexpr)
names(DEGresSig)
MRscore<-list()
for (i in seq(length(DEGresSig))) {
  MRscore[[i]]<-expr2MRscore(expr = myexpr[[i]],DEexpr = DEGresSig[[i]])
  names(MRscore)[i]<-names(myexpr)[i]
}
head(MRscore$Bacterial_infection$MRscore_value)
mergeMRscore<-bind_rows(MRscore[[1]]$MRscore_value,
                        MRscore[[2]]$MRscore_value,
                        MRscore[[3]]$MRscore_value,
                        MRscore[[4]]$MRscore_value,
                        MRscore[[5]]$MRscore_value,
                        MRscore[[6]]$MRscore_value)

head(mergeMRscore)
head(meta)
MRscore_meta<-merge(meta[,c(1,5,6)],mergeMRscore,by="sampleID")
head(MRscore_meta)
MRscore_meta$HeathyGroup<-MRscore_meta$Healthy0.Infection
MRscore_meta$HeathyGroup<-ifelse(MRscore_meta$Healthy0.Infection==0,"Healthy","Other")
MRscore_meta$BacGroup<-MRscore_meta$Type
MRscore_meta$BacGroup<-ifelse(MRscore_meta$Type!="Bacterial_infection","Other",MRscore_meta$BacGroup)
MRscore_meta$BacGroup<-ifelse(MRscore_meta$Type=="Healthy","Healthy",MRscore_meta$BacGroup)
MRscore_meta$VirGroup<-MRscore_meta$Type
MRscore_meta$VirGroup<-ifelse(MRscore_meta$Type!="Viral_infection","Other",MRscore_meta$VirGroup)
MRscore_meta$VirGroup<-ifelse(MRscore_meta$Type=="Healthy","Healthy",MRscore_meta$VirGroup)
MRscore_meta$TotalGroup<-MRscore_meta$BacGroup
MRscore_meta$TotalGroup<-ifelse(MRscore_meta$HeathyGroup=="Healthy","Healthy",MRscore_meta$TotalGroup)
MRscore_meta$TotalGroup<-ifelse(MRscore_meta$BacGroup=="Bacterial_infection","Bacterial_infection",MRscore_meta$TotalGroup)
MRscore_meta$TotalGroup<-ifelse(MRscore_meta$VirGroup=="Viral_infection","Viral_infection",MRscore_meta$TotalGroup)
head(MRscore_meta)
newdata<-MRscore_meta[,c(2,6:10)]
head(newdata)
#commparisons
library(ggsci)
library(ggpubr)
my_comparisons<-list(c("Healthy", "Bacterial_infection"),
                     c("Healthy", "Mixed_infection"),
                     c("Healthy", "Sepsis"),
                     c("Healthy", "Sepsis_surgical"),
                     c("Healthy", "Uninfected_inflammation"),
                     c("Healthy", "Viral_infection"))
my_comparisons1<-list(c("Healthy", "Other"))
my_comparisons2<-list(c("Healthy", "Bacterial_infection"),
                     c("Healthy", "Other"),
                     c("Bacterial_infection", "Other"))
my_comparisons3<-list(c("Healthy", "Viral_infection"),
                      c("Healthy", "Other"),
                      c("Viral_infection", "Other"))
my_comparisons4<-list(c("Healthy", "Viral_infection"),
                      c("Healthy", "Other"),
                      c("Viral_infection", "Other"),
                      c("Healthy", "Bacterial_infection"),
                      c("Healthy", "Other"),
                      c("Bacterial_infection", "Other"),
                      c("Viral_infection","Bacterial_infection"))

fun_to_plot <- function(data, group, variable,comparisons) {
  p <- ggboxplot(data, x=group, y=variable,fill = group, 
                #palette = c("#00AFBB", "#E7B800", "#FC4E07"), 
                add = "jitter")+
    stat_compare_means(comparisons = comparisons,label.y = c(0.8, 1,1.2,1.4,1.6,1.8,2))+
    scale_fill_aaas()+
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
          axis.title.x = element_blank())
  return(p)
}
col31<-c("#303841","#D72323","#377F5B","#375B7F","#F2FCFC","#f0027f",
         "#FAF8DE","#666666","#BDF1F6","#023782","#5e4fa2","#F1C40F",
         "#ff7f00","#cab2d6","#240041","#ffff99","#0E3BF0","#a65628",
         "#f781bf","#808FA6","#2EB872","#F0FFE1","#F33535","#011F4E",
         "#82B269","#D3C13E","#3F9DCD","#014E1F","#AFFFDF","#3D002E",
         "#3A554A")
barplot(rep(1,times=length(col31)),col=col31,border=cm.colors(length(col31)),
        axes=FALSE, main="cm.colors"); box()
colnames(newdata)
p<-fun_to_plot(newdata,group = "Type",variable = "MRscore",comparisons = my_comparisons)
p1<-fun_to_plot(newdata,group = "HeathyGroup",variable = "MRscore",comparisons = my_comparisons1)
p2<-fun_to_plot(newdata,group = "BacGroup",variable = "MRscore",comparisons = my_comparisons2)
p3<-fun_to_plot(newdata,group = "VirGroup",variable = "MRscore",comparisons = my_comparisons3)
p4<-fun_to_plot(newdata,group = "TotalGroup",variable = "MRscore",comparisons = my_comparisons4)
p+scale_fill_manual(values = c("#0E3BF0","#AFFFDF","#ff7f00","#375B7F","#F2FCFC","#D72323","#FAF8DE"))
p1+scale_fill_manual(values = c("#0E3BF0","#377F5B"))
p2+scale_fill_manual(values = c("#0E3BF0","#377F5B","#ff7f00"))
p3+scale_fill_manual(values = c("#0E3BF0","#D72323","#377F5B"))
p4+scale_fill_manual(values = c("#0E3BF0","#D72323","#377F5B","#ff7f00"))
