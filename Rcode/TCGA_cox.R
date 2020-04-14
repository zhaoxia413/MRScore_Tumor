library(survival)
library(survminer)
library(reshape2)
library(data.table)
library(stringr)
library(reshape2)
library(dplyr)
library(plyr)
library(survivalROC)
library(forestplot)
options(stringsAsFactors = F)
load(file = "../dataset/TCGA_results/TCGAMRegnesHeatmap.Rdata")
cli<-fread("../dataset/TCGA_data/TCGA_cli_new.csv")%>%as.data.frame()
cli[1:3,1:3]
MRscoreall[1:3,1:3]
MRscoreall$sampleID<-str_sub(MRscoreall$sampleID,1,12)
meta<-merge(MRscoreall,cli[,-2],by="sampleID")
colnames(exprall)<-str_sub(colnames(exprall),1,12)
exprall[1:3,1:3]
colnames(meta)

levels(factor(meta$Treatment_outcome))
meta$Treatment_outcome<-ifelse(meta$Treatment_outcome=="Complete Remission/Response","CR",
                               ifelse(meta$Treatment_outcome=="Partial Remission/Response" ,"PR",
                               ifelse(meta$Treatment_outcome=="Persistent Disease","PD",
                               ifelse(meta$Treatment_outcome=="Progressive Disease","PD",
                               ifelse(meta$Treatment_outcome=="Stable Disease","SD","NA")))))
 meta$Response<-ifelse(meta$Treatment_outcome=="PD","NR",
                      ifelse(meta$Treatment_outcome=="NA","NA","R"))
 levels(factor(meta$Treatment_outcome))
 levels(factor(meta$Response))
 #癌症零期（Stage 0）：检测到异常细胞，但异常细胞没有扩散到别处，这个时期也称为原位癌（carcinoma in situ）
 meta$Stage<-meta$AJCC_stage
 meta$Stage<-gsub("Stage III[A-C]","High_stage",meta$Stage)
 meta$Stage<-gsub("Stage IV[A-C]","High_stage",meta$Stage)
 meta$Stage<-gsub("Stage II[A-C]","Low_stage",meta$Stage)
 meta$Stage<-gsub("Stage I[A-C]","Low_stage",meta$Stage)
 meta$Stage<-gsub("Stage III","High_stage",meta$Stage)
 meta$Stage<-gsub("Stage IV","High_stage",meta$Stage)
 meta$Stage<-gsub("Stage II","Low_stage",meta$Stage)
 meta$Stage<-gsub("Stage I","Low_stage",meta$Stage)
 meta$Stage<-ifelse(meta$Stage=="High_stage"|meta$Stage=="Low_stage",meta$Stage,"NA")
 levels(factor(meta$Stage)) 
 levels(factor(meta$Gender))
 colnames(meta)
 data<-subset(meta,Stage!="NA")
 g2<-surv_cutpoint(
   data,
   time = "OStime",
   event = "OS",
   variables = c("Age")
 )
 summary(g2)
plot(g2)#Age==65
data$MRscore<-ifelse(data$MRscore>=median(data$MRscore),"High_MRscore","Low_MRscore")
data$Age<-ifelse(data$Age>=65,"High_age","Low_age")
coxfit<-coxph(Surv(OStime,OS) ~ MRscore+Stage+Age+Gender,data = data) 
p<-forest_model(coxfit,factor_separate_line=F, 
                 format_options = list(colour= "black", 
                                       shape = 12, 
                                       text_size = 4, 
                                       banded = T), 
                 theme = theme_forest())
 p
 covariates<-list("MRscore","Stage","Age","Gender")

 univ_formulas<-list()
 for (i in seq_len(length(covariates))) {
   univ_formulas[[i]] <- sapply(covariates[[i]],
                                function(x){
                                  ##print(x)
                                  as.formula(paste('Surv(OStime, OS)~', x))
                                })
 }
 univ_models<-list()
 univ_results<-list()
 res_single<-list()
 single_pick<-list()
 for (i in seq_along(univ_formulas)) {
     univ_models[[i]]<-coxph(univ_formulas[[i]][[1]], data = data)  
   names(univ_models)[i]<-covariates[i]
 }
 univ_results <- lapply(univ_models,
                               function(x){
                                 x <- summary(x)
                                 p.value =signif(x$waldtest["pvalue"], digits = 2)
                                 beta <- signif(x$coefficients[1], digits = 2)
                                 HR <- signif(x$coefficients[2], digits = 2)
                                 HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                                 HR.confint.upper <- signif(x$conf.int[,"upper .95"], 2)
                                 HR <- paste0(HR, " (",
                                              HR.confint.lower, "-", HR.confint.upper, ")")
                                 res <- c(beta, HR, p.value)
                                 names(res) <- c("coef", "HR (95% CI for HR)", "p.value")
                                 return(res)
                               })
 
 for (i in seq_along(univ_results)) {
   print(length(univ_results[[i]]))
 }
 
res_single <- as.data.frame(t(do.call(cbind, univ_results)))
table(res_single$p.value <= 0.05)
res_single <- res_single[order(res_single$p.value), ]
rownames(res_single) <- gsub('Surv.time..status....',"",rownames(res_single[[i]]))
single_pick <-rownames(res_single)
write.csv(res_single,"../dataset/TCGA_cox/KM_res_single_TCGA.csv")
fmla<- as.formula(paste0("Surv(OStime, OS) ~",paste0(single_pick,collapse = '+')))
cox<- coxph(fmla, data = data)
mutil_res<-summary(cox)
mutil_res_table<-as.data.frame(cbind(mutil_res$coefficients,mutil_res$conf.int[,-c(1,2)]))
rownames(mutil_res_table)<-gsub("low","",rownames(mutil_res$coefficients))
mutil_res_table<-subset(mutil_res_table,`Pr(>|z|)`<=0.05)
multi_pick<-rownames(mutil_res_table)
write.csv(mutil_res_table,"../dataset/TCGA_cox/mutil_res_table.csv")
final_res<-summary(cox)
cox=step(coxfit,direction = "both")
riskScore=predict(cox,type="risk",newdata=data)
risk=as.vector(ifelse(riskScore>1,"high","low"))
multiCOX_risk_result<-cbind(id=rownames(cbind(data[,8:9],riskScore,risk)),
                            cbind(data[,8:9],riskScore,risk))
write.csv(multiCOX_risk_result,"../dataset/TCGA_cox/mutation_model-multiCOX_risk.csv")
rt=read.csv("../dataset/TCGA_cox/mutation_model-multiCOX_risk.csv",check.names=F,row.names=1)
roc=survivalROC(Stime=rt$OStime, status=rt$OS, marker = rt$riskScore, 
                      predict.time =3, method="KM")
pdf(file=paste0("../dataset/TCGA_cox/mutation_model-multiCOX_riskROC.pdf"))
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
      xlab="False positive rate", ylab="True positive rate",
      main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
      lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
 abline(0,1)
 dev.off()
 diff=survdiff(Surv(OStime, OS) ~risk,data = rt)
 pValue=1-pchisq(diff$chisq,df=1)
 pValue=round(pValue,14)
 fit <- survfit(Surv(OStime/30, OS) ~ risk, data = rt)
 C_index<-data.frame(concordance=summary(cox)$concordance)
 pdf(file=paste0("../dataset/TCGA_cox/mutation_model-multiCOX_survival_curve.pdf"))
 plot(fit, lty = 2:3,col=c("red","blue"),xlab="time (months)",
      ylab="survival rate",
      main=paste("survival curve (p=", pValue[[i]] ,")",sep=""),mark.time=T)
 legend("topright", c("Low risk", "High risk"), lty = 2:3, col=c("blue","red"))
 dev.off()
 types<-levels(factor(meta$Types))
 df=subset(meta,Types==types[1])
 library(forestmodel)
 fun2forest<-function(data){
   print(i)
   coxfit<-coxph(Surv(OS.time,OS) ~ MRscore+Stage+Age+Gender,data = data) 
   p<-forest_model(coxfit,factor_separate_line=F, 
                   format_options = list(colour= "black", 
                                         shape = 12, 
                                         text_size = 4, 
                                         banded = T), 
                   theme = theme_forest())
   return(p)
 }
 p<-list()
 df<-list()
 data<-subset(meta,Stage!="NA")
 for (i in seq_along(types)) {
   tryCatch(
     {df[[i]]=subset(data,Types==types[i])
   df[[i]]$MRscore<-ifelse(df[[i]]$MRscore>=median(df[[i]]$MRscore),"High_MRscore","Low_MRscore")
   p[[i]]<-fun2forest(df[[i]])
   names(p)[i]<-types[i]
   pdf(file = paste0("../dataset/TCGA_cox/forest_",types[i],".pdf"),onefile = FALSE)
   print(p[[i]])
   dev.off()
   png(file = paste0("../dataset/TCGA_cox/forest_",types[i],".png"),width = 749,height = 364)
   print(p[[i]])
   dev.off()},
   error=function(e) NA)}
 
 "Surv(OS.time,OS) ~ MRscore+Stage+Age+Geender"
 fmla_BLCA=as.formula("Surv(OStime,OS) ~ MRscore+Stage+Age")
 fmla_BRCA=as.formula("Surv(OStime,OS) ~ MRscore+Stage+Age")
 fmla_HNSC=as.formula("Surv(OStime,OS) ~ MRscore+Stage+Age")
 fmla_LIHC=as.formula("Surv(OStime,OS) ~ MRscore+Stage")
 fmla_LUAD=as.formula("Surv(OStime,OS) ~ MRscore+Stage")
 
 fmlalist<-list(BLCA=fmla_BLCA,
            BRCA=fmla_BRCA,
            HNSC=fmla_HNSC,
            LIHC=fmla_LIHC,
            LUAD=fmla_LUAD)
 
 names(fmlalist)
 fmla=fmla_BLCA
 data<-subset(meta,Stage!="NA")%>%subset(.,Types=="BLCA")
 fmla2multifit<-function(fmla,data){
    cox<- coxph(fmla, data = data)
    mutil_res<-summary(cox)
    mutil_res_table<-as.data.frame(cbind(mutil_res$coefficients,mutil_res$conf.int[,-c(1,2)]))
    #mutil_res_table<-subset(mutil_res_table,`Pr(>|z|)`<=0.05)
    final_res<-summary(cox)
    riskScore=predict(cox,type="risk",newdata=data)
    risk=as.vector(ifelse(riskScore>1,"high","low"))
    multiCOX_risk_result<-cbind(id=rownames(cbind(data[,8:9],riskScore,risk)),
                             cbind(data[,8:9],riskScore,risk))
    rt=multiCOX_risk_result
    roc=survivalROC(Stime=rt$OStime, status=rt$OS, marker = rt$riskScore, 
                 predict.time =3, method="KM")
    par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
    rocplot<-plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
                  xlab="False positive rate", ylab="True positive rate",
                  main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
                  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
    abline(0,1)
    diff=survdiff(Surv(OStime, OS) ~risk,data = rt)
    pValue=1-pchisq(diff$chisq,df=1)
    pValue=round(pValue,14)
    fit <- survfit(Surv(OStime/30, OS) ~ risk, data = rt)
    C_index<-data.frame(concordance=summary(cox)$concordance)
    suvplot<-plot(fit, lty = 2:3,col=c("red","blue"),xlab="time (months)",
      ylab="survival rate",
      main=paste("survival curve (p=", pValue ,")",sep=""),mark.time=T)
    legend("topright", c("Low risk", "High risk"), lty = 2:3, col=c("blue","red"))
    return(list(multiCOX_risk_result=multiCOX_risk_result,
                mutil_res_table=mutil_res_table,
                rocplot=rocplot,
                suvplot=suvplot,
                C_index=C_index))
    }
 
multires<-list() 
data1<-list()
for (i in seq_along(fmlalist)) {
   print(i)
   data1[[i]]<-subset(meta,Stage!="NA")%>%subset(.,Types==names(fmlalist)[i])
  multires[[i]]<-fmla2multifit(fmlalist[[i]],data1[[i]])
  names(multires)[i]=names(fmlalist)[i]
  write.csv(multires[[i]]$multiCOX_risk_result,paste0("../dataset/TCGA_cox/mutation_model-multiCOX_risk_",names(fmlalist)[i],".csv"))
}
 
multires$BLCA$mutil_res_table
multires$BRCA$mutil_res_table
multires$HNSCmutil_res_table
multires$LIHC$mutil_res_table
multires$LUAD$mutil_res_table
multires$BLCA$suvplot
multires$BRCA$rocplot
multires$HNSC$rocplot
multires$LIHC$rocplot
multires$LUAD$rocplot
