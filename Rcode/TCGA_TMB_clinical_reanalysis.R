library(data.table)
library(survival)
library(survminer)
library(survey)
library(ggpubr)
library(ggsci)
library(ggthemes)
library(pheatmap)
tmb_micro<-fread("../dataset/TCGA_data/MRscore_TMB_patients7832.csv")%>%as.data.frame()
microScore<-fread("../dataset/TCGA_data/MRscore_TCGA_Patients9359.csv")%>%as.data.frame()
clinic<-fread("../dataset/TCGA_data/TCGA_clinical.csv")%>%as.data.frame()
colnames(clinic)[1]<-"Patient_ID"
colnames(microScore)[1]<-"Patient_ID"
microScore$Patient_ID<-str_sub(microScore$Patient_ID,1,12)
microScore$Patient_ID[1:3]
data<-merge(tmb_micro,clinic,by="Patient_ID")
head(microScore)
clinic[1:3,1:3]
data1<-merge(microScore,clinic,by="Patient_ID")
##re plot the survival curve of tcga
datalist<-split.data.frame(data1,f=data1$Types,drop = F)
surv_cut_OS<-list()
fit<-list()
p<-list()
###OS
seq_along(datalist)
try(for (i in c(1:20,22,23,25:28)) {
  print(i)
  surv_cut_OS[[i]] <- surv_cutpoint(
    datalist[[i]],
    time = "OS.time",
    event = "OS",
    variables = c("MRscore")
  )
  summary(surv_cut_OS[[i]])
  plot(surv_cut_OS[[i]])
  datalist[[i]]$MRscore_group<-sample(c("High","Low"),nrow(datalist[[i]]),replace = T)
  datalist[[i]]$MRscore_group<-ifelse(datalist[[i]]$MRscore>=summary(surv_cut_OS[[i]])[,1],"High","Low")
  fit[[i]]<-survfit(Surv(OS.time,OS) ~ MRscore_group,
                    data = datalist[[i]])
  p[[i]]<-ggsurvplot(fit[[i]],
                     data=datalist[[i]],
                     risk.table = TRUE,
                     pval = TRUE,
                     ggtheme = theme_survminer(),
                     title = paste0("OS_",names(datalist)[i]),
                     palette = c("blue","red"),
                     #facet.by = "Types",
                     legend.title="MRscore_",
                     risk.table.col = "strata",
                     surv.median.line = "hv",
                     risk.table.y.text.col = T,
                     risk.table.y.text = FALSE )
  pdf(file = paste0("../dataset/TCGA_results/survical_MRscore/",names(datalist)[i],"_OS.pdf"),width = 5,height = 6,onefile = FALSE)
  print(p[[i]])
  dev.off()
  png(file = paste0("../dataset/TCGA_results/survical_MRscore/",names(datalist)[i],"_OS.png"))
  print(p[[i]])
  dev.off()
},silent = T)

##PFS
for (i in c(1:13,15:20,22,23,25:28)) {
  print(i)
  surv_cut_OS[[i]] <- surv_cutpoint(
    datalist[[i]],
    time = "PFS.time",
    event = "PFS",
    variables = c("MRscore")
  )
  summary(surv_cut_OS[[i]])
  plot(surv_cut_OS[[i]])
  datalist[[i]]$MRscore_group<-sample(c("High","Low"),nrow(datalist[[i]]),replace = T)
  datalist[[i]]$MRscore_group<-ifelse(datalist[[i]]$MRscore>=summary(surv_cut_OS[[i]])[,1],"High","Low")
  fit[[i]]<-survfit(Surv(PFS.time,PFS) ~ MRscore_group,
                    data = datalist[[i]])
  p[[i]]<-ggsurvplot(fit[[i]],
                     data=datalist[[i]],
                     risk.table = TRUE,
                     pval = TRUE,
                     ggtheme = theme_survminer(),
                     title = paste0("PFS_",names(datalist)[i]),
                     palette = c("blue","red"),
                     #facet.by = "Types",
                     legend.title="MRscore_",
                     risk.table.col = "strata",
                     surv.median.line = "hv",
                     risk.table.y.text.col = T,
                     risk.table.y.text = FALSE )
  pdf(file = paste0("../dataset/TCGA_results/survical_MRscore/",names(datalist)[i],"_PFS.pdf"),width = 5,height = 6,onefile = FALSE)
  print(p[[i]])
  dev.off()
  png(file = paste0("../dataset/TCGA_results/survical_MRscore/",names(datalist)[i],"_PFS.png"))
  print(p[[i]])
  dev.off()
}
#TMB frequence
library(RColorBrewer)
library(ggridges)

colnames(data)[5]<-"TMB"
subdata<-subset(data,TMB<=8)#data 7827samples subdata 7090
subdata$Types<-subdata$Types[order(subdata$Types,decreasing = F)]
Freq<-ggplot(subdata, 
             aes(TMB,Types, fill = ..density..)) + 
  geom_density_ridges_gradient(aes(height = ..density..),scale = 1,size = 0.3)+#sacle���÷�???
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11,'BrBG')))(32))+
  labs(x="TMB")+
  theme_few(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5))
#facet_wrap(~Immune_Subtype,scales = "free_x")
Freq
Freq1<-ggplot(subdata, 
             aes(MRscore,Types, fill = ..density..)) + 
  geom_density_ridges_gradient(aes(height = ..density..),scale = 1,size = 0.3)+#sacle���÷�???
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(9,'YlGnBu')))(32))+
  labs(x="MRscore")+
  theme_few(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5))
#facet_wrap(~Immune_Subtype,scales = "free_x")
Freq1
##TMB_MRscore_OS_cutoff
datalist1<-split.data.frame(subdata,f=subdata$Types,drop = F)
fun_to_plot <- function(data, group, variable) {
  p <- ggboxplot(data, x=group, y=variable,fill = group, 
                 palette = c("#00AFBB", "#E7B800", "#FC4E07"), 
                 add = "jitter", shape=group)+
    stat_compare_means(comparisons = my_comparisons)
  return(p)
}
surv_cut_OS<-list()
p<-list()
###OS
seq_along(datalist1)
my_comparisons<-list(c("High","Low"))
for (i in seq_along(datalist1)) {
  print(i)
  surv_cut_OS[[i]] <- surv_cutpoint(
    datalist1[[i]],
    time = "OS.time",
    event = "OS",
    variables = c("MRscore")
  )
  summary(surv_cut_OS[[i]])
  plot(surv_cut_OS[[i]])
  datalist1[[i]]$MRscore_level<-sample(c("High","Low"),nrow(datalist1[[i]]),replace = T)
  datalist1[[i]]$MRscore_level<-ifelse(datalist1[[i]]$MRscore>=summary(surv_cut_OS[[i]])[,1],"High","Low")
  datalist1[[i]]$MRscore_level<-factor( datalist1[[i]]$MRscore_level,levels = c("High","Low"))
  p[[i]]<-fun_to_plot(datalist1[[i]],"MRscore_level","TMB")
  pdf(file = paste0("../dataset/TCGA_results/survical_MRscore/",names(datalist1)[i],"_OScutMRscore_TMB.pdf"),width = 5,height = 6,onefile = FALSE)
  print(p[[i]])
  dev.off()
  png(file = paste0("../dataset/TCGA_results/survical_MRscore/",names(datalist1)[i],"_OScutMRscore_TMB.png"))
  print(p[[i]])
  dev.off()
  }
##PFS_cutoff MRscore level
for (i in c(1:11,13:22)) {
  print(i)
  surv_cut_OS[[i]] <- surv_cutpoint(
    datalist1[[i]],
    time = "PFS.time",
    event = "PFS",
    variables = c("MRscore")
  )
  summary(surv_cut_OS[[i]])
  plot(surv_cut_OS[[i]])
  datalist1[[i]]$MRscore_level<-sample(c("High","Low"),nrow(datalist1[[i]]),replace = T)
  datalist1[[i]]$MRscore_level<-ifelse(datalist1[[i]]$MRscore>=summary(surv_cut_OS[[i]])[,1],"High","Low")
  datalist1[[i]]$MRscore_level<-factor( datalist1[[i]]$MRscore_level,levels = c("High","Low"))
  p[[i]]<-fun_to_plot(datalist1[[i]],"MRscore_level","TMB")
  pdf(file = paste0("../dataset/TCGA_results/survical_MRscore/",names(datalist1)[i],"_PFScutMRscore_TMB.pdf"),width = 5,height = 6,onefile = FALSE)
  print(p[[i]])
  dev.off()
  png(file = paste0("../dataset/TCGA_results/survical_MRscore/",names(datalist1)[i],"_PFScutMRscore_TMB.png"))
  print(p[[i]])
  dev.off()
}
plotdata<-bind_rows(datalist1)
p <-ggplot(plotdata,aes(Types,TMB,fill = MRscore_level)) +
  geom_boxplot()+
  scale_fill_manual(values = c("#00AFBB", "#E7B800"))+
  theme_bw()+
  coord_flip()
p  
p <-ggplot(plotdata,aes(Types,TMB,fill = MRscore_level)) +
  geom_boxplot()+
  scale_fill_d3()+
  theme_bw()+
  coord_flip()
p  



