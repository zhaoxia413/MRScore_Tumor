library(data.table)
library(ggsci)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
options(stringsAsFactors = F)
df<-fread("../dataset/TCGA_data/MRgenes_CancerTYpes_DEGtag.csv")
up<-subset(df,Regulation=="Up")
down<-subset(df,Regulation=="Down")
genelist_up<-split(up,f=up$Types)
genelist_down<-split(down,f=down$Types)
keytypes(org.Hs.eg.db)
upgenes_ID<-list()
downgene_ID<-list()
upgenes<-list()
downgenes<-list()
ekkUp<-list()
ekkDown<-list()
for (i in seq(length(genelist_up))) {
  upgenes_ID[[i]]=select(org.Hs.eg.db,keys = genelist_up[[i]]$Gene,column = "ENTREZID", keytype = "SYMBOL",multiVals = "first")
  downgene_ID[[i]]=select(org.Hs.eg.db,keys = genelist_down[[i]]$Gene,column = "ENTREZID", keytype = "SYMBOL",multiVals = "first")
  upgenes[[i]]<-upgenes_ID[[i]][,2]
  downgenes[[i]]<-downgene_ID[[i]][,2]
  #organism:'species' should be one of organisms listed in 'http://www.genome.jp/kegg/catalog/org_list.html'
  ekkUp[[i]] <- enrichKEGG(gene = upgenes[[i]],organism = 'hsa',pvalueCutoff = 0.05)
  ekkDown[[i]] <- enrichKEGG(gene = downgenes[[i]],organism = 'hsa',pvalueCutoff = 0.05)
write.csv(as.matrix(ekkUp[[i]]@result),
          file=paste0("../dataset/TCGA_results/KEGGenrich/",names(genelist_up)[i],"_KEGG_up.txt"),
          row.names = F,sep = "\t",quote = F)
write.table(as.matrix(ekkDown[[i]]@result), 
            file=paste0("../dataset/TCGA_results/KEGGenrich/",names(genelist_down)[i],"_KEGG_Down.txt"),
            row.names = F,sep = "\t",quote = F)
}
files<-list.files("../dataset/TCGA_results/KEGGenrich/")
upfiles<-files[grep("up.txt$",files)]
downfiles<-files[grep("Down.txt$",files)]
uppath<-paste0("../dataset/TCGA_results/KEGGenrich/",upfiles)
downpath<-paste0("../dataset/TCGA_results/KEGGenrich/",downfiles)
kegg_up<-list()
kegg_down<-list()
for (i in 1:28) {
  kegg_up[[i]]<-fread(uppath[i],
                      select = c("Description","GeneRatio","BgRatio","qvalue","Count"))%>%
    subset(.,qvalue<=0.05)
  kegg_up[[i]]$GeneRatio<-as.numeric(gsub("/.*","",kegg_up[[i]]$GeneRatio))
  kegg_up[[i]]$BgRatio<-as.numeric(gsub("/.*","",kegg_up[[i]]$BgRatio))
  kegg_up[[i]]$Ratio=kegg_up[[i]]$GeneRatio/kegg_up[[i]]$BgRatio
  kegg_up[[i]]$negLogqvalue=-(log(kegg_up[[i]]$qvalue))
  kegg_up[[i]]<-kegg_up[[i]][,-c(2:4)]
  kegg_up[[i]]$regulation<-rep("Up",nrow(kegg_up[[i]]))
  names(kegg_up)[i]<-gsub("_KEGG.*","",upfiles[i])
  kegg_up[[i]]$CancerTypes<-rep(names(kegg_up)[i],nrow(kegg_up[[i]]))
  kegg_up[[i]]$enrichScore=log(kegg_up[[i]]$Count+kegg_up[[i]]$Ratio+kegg_up[[i]]$negLogqvalue)
  kegg_up[[i]]<-top_n( kegg_up[[i]],5,enrichScore)
  kegg_down[[i]]<-fread(downpath[i],
                      select = c("Description","GeneRatio","BgRatio","qvalue","Count"))%>%
    subset(.,qvalue<=0.05)
  kegg_down[[i]]$GeneRatio<-as.numeric(gsub("/.*","",kegg_down[[i]]$GeneRatio))
  kegg_down[[i]]$BgRatio<-as.numeric(gsub("/.*","",kegg_down[[i]]$BgRatio))
  kegg_down[[i]]$Gene<-as.numeric(gsub("*./","",kegg_down[[i]]$GeneRatio))
  kegg_down[[i]]$Bg<-as.numeric(gsub("*./","",kegg_down[[i]]$BgRatio))
  kegg_down[[i]]$Ratio=kegg_down[[i]]$GeneRatio/kegg_down[[i]]$BgRatio
  kegg_down[[i]]$negLogqvalue=-(log(kegg_down[[i]]$qvalue))
  kegg_down[[i]]<-kegg_down[[i]][,-c(2:4)]
  kegg_down[[i]]$regulation<-rep("Down",nrow(kegg_down[[i]]))
  names(kegg_down)[i]<-gsub("_KEGG.*","",downfiles[i])
  kegg_down[[i]]$CancerTypes<-rep(names(kegg_down)[i],nrow(kegg_down[[i]]))
  kegg_down[[i]]$enrichScore=log(kegg_down[[i]]$Count+kegg_down[[i]]$Ratio+kegg_down[[i]]$negLogqvalue)
  kegg_down[[i]]<-top_n( kegg_down[[i]],5,enrichScore)
}
sapply(kegg_up, nrow)
sapply(kegg_down, nrow)
kegg_up<-kegg_up[-which(sapply(kegg_up, nrow)==0)]
kegg_down<-kegg_down[-which(sapply(kegg_down, nrow)==0)]
sapply(kegg_up, nrow)
sapply(kegg_down, nrow)
keggall<-bind_rows(kegg_up,kegg_down)
write.table(keggall,'../dataset/TCGA_results/KEGGenrich/allKEGG.txt',row.names = F,sep = "\t",quote = F)

library(ggpubr)
library(ggthemes)
keggall<-fread("../dataset/TCGA_results/KEGGenrich/allKEGG.txt",header = T)
head(keggall)
kegg_summary<-keggall%>%group_by(Description,regulation)%>%summarise(counts=n())
head(kegg_summary)
kegg_summary<-subset(kegg_summary,abs(counts)>=5)
kegg_summary$Description<-reorder(kegg_summary$Description,kegg_summary$counts)
pk<-ggbarplot(kegg_summary, x="Description", y="counts", 
              fill = "regulation", color = "white", 
              width = 1,palette = "npg", 
               sort.by.groups = F, 
                          x.text.angle=90, ylab = "GeneNumber", 
                          xlab = "KEGG_term", 
              legend.title="Class", 
              rotate=TRUE, 
              ggtheme = theme_tufte(base_size = 12))
pk
keyTypes<-c("ACC","KIRC","KIRP","PRAD","STAD","THCA","BLCA","LUSC")
keggall1<-keggall[which(keggall$CancerTypes%in%keyTypes),]
keggall2<-keggall[-which(keggall$CancerTypes%in%keyTypes),]
keggall1$Description<-reorder(keggall1$Description,keggall1$enrichScore)
keggall2$Description<-reorder(keggall2$Description,keggall2$enrichScore)
p<-ggplot(keggall1,aes(enrichScore,Description))
p<-ggplot(keggall2,aes(enrichScore,Description))
pp<-p+geom_point(aes(size=Count,color=negLogqvalue,shape=regulation))
p1<-pp+ scale_color_gradient(low = "blue",high = "red")+
  theme(axis.text = element_text(colour = "black",size = 12),
        legend.text = element_text(colour = "black",size = 12),
        panel.grid.major = element_line(colour = "gray88",linetype = "dashed"),
        panel.background=element_blank(),
        axis.ticks = element_line(size = 1, colour = "black"),
        axis.line = element_line(size = 1, colour = "black"),
        legend.position = "right",
        panel.border = element_rect(linetype = "solid", fill = NA,size = 2),
        axis.title = element_text(size = 12)
        ,title = element_text(size = 12))+
  facet_wrap(~CancerTypes,scales = "free",ncol = 4)
p1+labs(x="enrichScore")+
  theme(axis.text = element_text(colour = "black",size = 10),
        legend.text = element_text(colour = "black",size = 10),
        panel.grid.major = element_line(colour = "gray88",linetype = "dashed"),
        panel.background=element_blank(),
        axis.ticks = element_line(size = 1, colour = "black"),
        axis.line = element_line(size = 1, colour = "black"),
        legend.position = "right",
        panel.border = element_rect(linetype = "solid", fill = NA,size = 2),
        axis.title = element_text(size = 12)
        ,title = element_text(size = 12))

##enrichGO

egoUp1<-list()
egoUp2<-list()
egoUp3<-list()
egoDown1<-list()
egoDown2<-list()
egoDown3<-list()
rUp1<-list()
rUp2<-list()
rUp3<-list()
rDown1<-list()
rDown2<-list()
rDown3<-list()
egoUp<-list()
egoDown<-list()
for (i in seq(length(genelist_up))) {
egoUp1[[i]] <- enrichGO(gene=genelist_up[[i]]$Gene,org.Hs.eg.db,ont="BP",pvalueCutoff=0.05,readable=TRUE)
egoDown1[[i]] <- enrichGO(gene=genelist_down[[i]]$Gene,org.Hs.eg.db,ont="BP",pvalueCutoff=0.05,readable=TRUE)
rUp1[[i]]<-rep("BP",length(row_number(egoUp1[[i]]@result$Description)))
rUp1[[i]]<-as.data.frame(rUp1[[i]])
colnames(rUp1[[i]])<-("Class")
egoUp1[[i]]<-bind_cols(egoUp1[[i]]@result,rUp1[[i]])

rDown1[[i]]<-rep("BP",length(row_number(egoDown1[[i]]@result$Description)))
rDown1[[i]]<-as.data.frame(rDown1[[i]])
colnames(rDown1[[i]])<-("Class")
egoDown1[[i]]<-bind_cols(egoDown1[[i]]@result,rDown1[[i]])

egoUp2[[i]] <- enrichGO(gene=geneUp,org.Hs.eg.db,ont="MF",pvalueCutoff=0.05,readable=TRUE)
egoDown2[[i]]<- enrichGO(gene=geneDown,org.Hs.eg.db,ont="MF",pvalueCutoff=0.05,readable=TRUE)

rUp2[[i]]<-rep("MF",length(row_number(egoUp2[[i]]@result$Description)))
rUp2[[i]]<-as.data.frame(rUp2[[i]])
colnames(rUp2[[i]])<-("Class")
egoUp2[[i]]<-bind_cols(egoUp2[[i]]@result,rUp2[[i]])

rDown2[[i]]<-rep("MF",length(row_number(egoDown2[[i]]@result$Description)))
rDown2[[i]]<-as.data.frame(rDown2[[i]])
colnames(rDown2[[i]])<-("Class")
egoDown2[[i]]<-bind_cols(egoDown2[[i]]@result,rDown2[[i]])

egoUp3[[i]] <- enrichGO(gene=geneUp,org.Hs.eg.db,ont="CC",pvalueCutoff=0.05,readable=TRUE)
egoDown3[[i]]<- enrichGO(gene=geneDown,org.Hs.eg.db,ont="CC",pvalueCutoff=0.05,readable=TRUE)

rUp3[[i]]<-rep("CC",length(row_number(egoUp3[[i]]@result$Description)))
rUp3[[i]]<-as.data.frame(rUp3[[i]])
colnames(rUp3[[i]])<-("Class")
egoUp3[[i]]<-bind_cols(egoUp3[[i]]@result,rUp3[[i]])

rDown3[[i]]<-rep("CC",length(row_number(egoDown3[[i]]@result$Description)))
rDown3[[i]]<-as.data.frame(rDown3[[i]])
colnames(rDown3[[i]])<-("Class")
egoDown3[[i]]<-bind_cols(egoDown3[[i]]@result,rDown3[[i]])

egoUp[[i]]<-bind_rows(egoUp1[[i]],egoUp2[[i]],egoUp3[[i]])
egoDown[[i]]<-bind_rows(egoDown1[[i]],egoDown2[[i]],egoDown3[[i]])
write.table(egoUp[[i]],
          file=paste0("../dataset/TCGA_results/GOenrich/",names(genelist_up)[i],"_GO_up.txt"),
          row.names = F,sep = "\t",quote = F)
write.table(egoDown[[i]], 
            file=paste0("../dataset/TCGA_results/GOenrich/",names(genelist_down)[i],"_GO_Down.txt"),
            row.names = F,sep = "\t",quote = F)
}
keypathway<-fread("../dataset/TCGA_data/Nod_TLR_genes.csv")
DE<-fread("../dataset/TCGA_data/MRgenes_CancerTYpes_DEGtag.csv")
keypathwaygeneExpr<-merge(keypathway,DE,by="Gene")
head(keypathwaygeneExpr)
voldata<-keypathwaygeneExpr[,c(1,4,5,6,7)]
voldata$regulation<-factor(ifelse(abs(voldata$Log2FC) >= 5, ifelse(voldata$Log2FC>= 5 ,'Up','Down'),'FC (1-5)'),levels=c('Up','Down','FC (1-5)'))
p<-ggplot(voldata,aes(Log2FC,-log10(adjp),
                                colour=regulation))+
  ylab("-log10(adjp)")+
  geom_point(aes(shape=KEGG_term,size=abs(Log2FC)),alpha=0.6)+
  geom_vline(xintercept=c(-5,5),lty=4,col="black",lwd=0.8) +
  theme_classic(base_size = 16)+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+
  geom_text_repel(
    data = subset(voldata,abs(voldata$Log2FC)>=5),
    aes(label = Gene),
    color="black",
    size = 4,
    box.padding = unit(1,"lines"),
    point.padding = unit(1,"lines"),
    segment.color = "black",
    show.legend = F
  )
p
p1<-ggplot(voldata,aes(Log2FC,-log10(adjp),
                      colour=factor(Log2FC>0)))+
  ylab("-log10(adjp)")+
  geom_point(aes(shape=KEGG_term,size=abs(Log2FC)),alpha=0.6)+
  theme_classic(base_size = 16)+
  scale_color_aaas()+
  geom_text_repel(
    data = subset(voldata,abs(voldata$Log2FC)>=5),
    aes(label = Gene),
    color="black",
    size = 4,
    box.padding = unit(1,"lines"),
    point.padding = unit(1,"lines"),
    segment.color = "black",
    show.legend = F
  )
p1+facet_wrap(~Types,scales = "free")

voldata<-split(voldata,voldata$Types)
library(ggrepel)
library(grid)
library(gridExtra)
p2<-list()
for (i in seq(length(voldata))) {
  print(i)
p2[[i]]<-ggplot(voldata[[i]],aes(Log2FC,-log10(adjp),
                    colour=factor(Log2FC<0)))+
  ylab("-log10(adjp)")+
  geom_point(aes(shape=KEGG_term,size=abs(Log2FC)),alpha=0.6)+
  geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8)+
  theme_classic(base_size = 16)+
  scale_color_aaas()+
  ggtitle(label = names(voldata)[i])+
  geom_text_repel(
    data = top_n(voldata[[i]],5,abs(voldata[[i]]$Log2FC)),
    aes(label = Gene),
    color="black",
    size = 4,
    box.padding = unit(1,"lines"),
    point.padding = unit(1,"lines"),
    segment.color = "black",
    show.legend = F
  )
  grid.newpage()
  pdf(paste0("../dataset/TCGA_results/KEGGenrich/kegg_keygene",names(voldata)[i],".pdf"),width = 7,height = 5)
  grid.draw(p2[[i]])
  dev.off()
  grid.newpage()
  png(paste0("../dataset/TCGA_results/KEGGenrich/kegg_keygene",names(voldata)[i],".png"))
  grid.draw(p2[[i]])
  dev.off()
}

library(reshape2)
library(ggridges)
library(RColorBrewer)
head(keypathwaygeneExpr)
df<-subset(keypathwaygeneExpr,abs(Log2FC)>=5)
DE_TMP<-fread("../dataset/TCGA_data/microSignature_DEG_TMP.csv")
keygeneTMP<-merge(df[,c(1,4)],DE_TMP,by="Gene")
keygene<-levels(factor(df$Gene))
data<-keypathwaygeneExpr[,c(1,4,5,7)][which(Gene%in%keygene),]
head(data)
Freq<-ggplot(data, 
              aes(Log2FC,Gene, fill = ..density..)) + 
  geom_density_ridges_gradient(aes(height = ..density..),scale = 1,size = 0.3)+
  theme_minimal(base_size = 12)+
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11,'Spectral')))(32))+
  labs(x="FoldChange")+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~KEGG_term,scales = "free")
Freq+theme_cleveland()
library(pheatmap)
keygeneTMP<-keygeneTMP[!which(duplicated(keygeneTMP$Gene))]
keygeneTMP<-keygeneTMP[,-"V1"]
mat<-data.frame(row.names = keygeneTMP$Gene,keygeneTMP[,-c(1:2)])
sampleinfo<-read.csv("../dataset/TCGA_data/sampleInfo.csv",row.names = 1,header = T)
rownames(sampleinfo)<-gsub("-",".",rownames(sampleinfo))
setdiff(rownames(sampleinfo),colnames(mat))
pheatmap(scale(mat),show_colnames = F,annotation_col = sampleinfo,
         cluster_cols = F)
keyTypes<-c("ACC","KIRC","KIRP","PRAD","STAD","THCA","BLCA","LUSC")
sample1<-subset(sampleinfo,Types%in%keyTypes)
mat1<-mat[,which(colnames(mat)%in%rownames(sample1))]
pheatmap(scale(mat1),show_colnames = F,annotation_col = sample1,
         cluster_cols = T)
#ConsensusClusterPlusresults by keygene<-levels(factor(df$Gene))
library(ConsensusClusterPlus)
dir.create("../dataset/TCGA_data/ConsensusClusterPlusresults")
title="../dataset/TCGA_data/ConsensusClusterPlusresults"
d<-log(mat+1)
mads=apply(d,1,mad)
d = sweep(d,1, apply(d,1,median,na.rm=T))
results<-ConsensusClusterPlus(d,maxK=8,reps=50,pItem=0.8,pFeature=1,
                              title=title,clusterAlg="hc",distance="spearman",seed=1262118388.71279,plot="png")
cluster<-results[[3]]["consensusClass"]
write.csv(cluster,"../dataset/TCGA_data/ConsensusCluster.csv")
allDEfiles<-list.files("../../TCGA_DEGs/DEGresults/")
path<-paste0("../../TCGA_DEGs/DEGresults/",allDEfiles)
cancers<-str_sub(allDEfiles,15,18)
cancers[c(1,8,14,18,28)]<-c("ACC","GBM","LGG","OV","UCS")
allDE<-list()
MR_FC<-list()
keypathway<-fread("../dataset/TCGA_data/Nod_TLR_genes.csv",
                  select = c("Gene","KEGG_term"))
for (i in seq(length(allDEfiles))) {
  allDE[[i]]<-fread(path[i],header = T)
  colnames(allDE[[i]])[c(1,5)]<-c("Gene",cancers[i])
  allDE[[i]]<-allDE[[i]][,c(1,5)]
  MR_FC[[i]]<-left_join(keypathway,allDE[[i]])
}
MR_FC[[24]]<-MR_FC[[24]][-203,]
sapply(MR_FC, nrow)
MRgenes_FC<-bind_cols(MR_FC)
MRgenes_FC<-MRgenes_FC[,-grep("term[1-9]",colnames(MRgenes_FC))]
MRgenes_FC<-MRgenes_FC[,-grep("Gene[1-9]",colnames(MRgenes_FC))]
MRgenes_FC[is.na(MRgenes_FC)]=0
head(MRgenes_FC)
MRgenes_FC<-MRgenes_FC[!duplicated(MRgenes_FC$Gene),]
keypathway<-fread("../dataset/TCGA_data/Nod_TLR_genes.csv")
DE<-fread("../dataset/TCGA_data/MRgenes_CancerTYpes_DEGtag.csv")
keypathwaygeneExpr<-merge(keypathway,DE,by="Gene")
DE_TMP<-fread("../dataset/TCGA_data/microSignature_DEG_TMP.csv")
df<-subset(keypathwaygeneExpr,abs(Log2FC)>=5)
keygeneTMP<-merge(df[,c(1,4)],DE_TMP,by="Gene")
keygene<-levels(factor(df$Gene))
MRgenes_FC<-MRgenes_FC[MRgenes_FC$Gene%in%keygene,]
mat<-data.frame(row.names = MRgenes_FC$Gene,MRgenes_FC[,-c(1:2)])
library(viridis)
pheatmap(scale(mat,center = F),border_color = "NA",fontsize = 12)
