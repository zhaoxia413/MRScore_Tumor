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
<<<<<<< HEAD
 ##glmnet continus variables cox
 library(glmnet)
 library(survival)
 d <- subset(meta,select = c(OStime,OS,MRscore,Age,Stage))
 d<-subset(d,Age!="NA")%>%subset(.,OStime!="NA")
 d$OStime=round(d$OStime/30,2)%>%as.numeric()
 y <- Surv(d$OStime, d$OS)
 x <- model.matrix(y ~ MRscore + Age + Stage, d)
 fit <- glmnet(x, y, alpha=1,family="mgaussian")
 plot(fit, label=T)
 cv.fit <- cv.glmnet(x, y, family="cox", alpha=1)
 cv_fit <- cv.glmnet(x, y, nfolds = 5) #just with numeric variables
 cv_fit_2 <- cv.glmnet(cbind(x ,x_train), y, nfolds = 5) #both factor and numeric variables
 par(mfrow = c(2,1))
 plot(cv_fit)
 plot(cv_fit_2)
 plot(cv.fit)
=======

>>>>>>> 6255cae2aadc6bf9f411b1ac81dc243fbe02c9cf
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
 colnames(meta)[c(9,11,13)]<-c("OStime","DSStime","PFItime")
 meta$Stage_ct<-meta$AJCC_stage
 meta$Stage_ct<-gsub("Stage III[A-C]","Stage III",meta$Stage_ct)
 meta$Stage_ct<-gsub("Stage IV[A-C]","Stage IV",meta$Stage_ct)
 meta$Stage_ct<-gsub("Stage II[A-C]","Stage II",meta$Stage_ct)
 meta$Stage_ct<-gsub("Stage I[A-C]","Stage I",meta$Stage_ct)
 meta$Stage_ct<-gsub("Stage III","Stage III",meta$Stage_ct)
 meta$Stage_ct<-gsub("Stage IV","Stage IV",meta$Stage_ct)
 meta$Stage_ct<-gsub("Stage II","Stage II",meta$Stage_ct)
 meta$Stage_ct<-gsub("Stage I","Stage I",meta$Stage_ct)
 meta$Stage_ct<-ifelse(meta$Stage_ct%in%c("Stage I","Stage II","Stage III","Stage IV"),meta$Stage_ct,"NA")
 levels(factor(meta$Stage)) 
 levels(factor(meta$Gender))
 levels(factor(meta$Stage_ct))
 ##reponse

 fun_to_plot <- function(data, group, variable,facet,comparisons) {
   col16<-c("#b3e2cd","#fdcdac","#cbd5e8","#1b9e77","#d95f02","#7570b3",
            "#e7298a","#66a61e","#e6ab02","#a6761d","#666666",
            "#f4cae4","#e6f5c9","#fff2ae","#f1e2cc","#cccccc")
   p <- ggboxplot(data, x=group, y=variable,fill = group, facet.by = facet,font.label = list(size=6),
                  palette = col16[sample(c(1:13),10)])+ #palette = col16[sample(c(1:10),10)]
     stat_compare_means(comparisons = comparisons)+
     geom_jitter(col=rgb(0, 0, 00, 0.2),size=1)
   return(p)
 }
 my_comparisons<-list(c("R","NR"))
 my_comparisons<-list(c("MALE","FEMALE"))
 my_comparisons<-list(c("High","Low"))
 my_comparisons<-list(c("Stage I","Stage II"),
                      c("Stage I","Stage III"),
                      c("Stage I","Stage IV"),
                      c("Stage II","Stage III"),
                      c("Stage II","Stage IV"),
                      c("Stage IV","Stage III"))
 
 data1=subset(meta,Response!="NA")
 data2=subset(meta,Response!="NA")%>%subset(.,Stage_ct!="NA")
 data2$Stage_ct<-factor(data2$Stage_ct,levels = c("Stage I","Stage II","Stage III","Stage IV"))
 data4=subset(meta,Gender!="NA")
 data5=subset(meta,Age!="NA")%>%subset(.,Stage_ct!="NA")
 fun_to_plot(data1,"Response","MRscore","Types",my_comparisons)
 fun_to_plot(data2,"Response","MRscore","Stage_ct",my_comparisons)
 fun_to_plot(data2,"Response","MRscore","Stage",my_comparisons)
 fun_to_plot(data1,"Response","MRscore",NULL,my_comparisons)
 fun_to_plot(data2,"Gender","MRscore","Stage_ct",my_comparisons)+theme(axis.text.x = element_text(angle = 90))
 fun_to_plot(data1,"Gender","MRscore",NULL,my_comparisons)
 fun_to_plot(data2,"Stage_ct","MRscore",NULL,my_comparisons)+theme(axis.text.x = element_text(angle = 90))
 fun_to_plot(data1,"Gender","MRscore",NULL,my_comparisons)
 data2$MRscore_group<-ifelse(data2$MRscore>=mean(meta$MRscore),"High","Low")
 data5$Age_60<-ifelse(data5$Age>=60,"High","Low")
 fun_to_plot(data5,"Age_60","MRscore","Stage_ct",my_comparisons)
 fun_to_plot(data5,"Age_60","MRscore",NULL,my_comparisons)
 
 fit<-survfit(Surv(OStime/365,OS) ~ MRscore_group,ctype = 1,
              conf.type="log-log",
              data = data2)
 surv_median(fit)
 ggsurvplot(fit = fit,facet.by = "Stage_ct",
            data = data2,
            pval.method = T,
            risk.table = TRUE,
            pval = TRUE,
            ggtheme = theme_survminer(),
            palette = c("red","blue"),
            #facet.by = "Types",
            legend.title="MRscore(mean)",
            risk.table.col = "strata",
            surv.median.line = "hv",
            risk.table.y.text.col = T,
            risk.table.y.text = T)
 fit<-survfit(Surv(PFItime/365,PFI) ~ MRscore_group,ctype = 1,
              conf.type="log-log",
              data = data2)
 surv_median(fit)
 ggsurvplot(fit = fit,facet.by = "Stage_ct",
            data = data2,
            pval.method = T,
            risk.table = TRUE,
            pval = TRUE,
            ggtheme = theme_survminer(),
            palette = c("#08d9d6","#252a34"),
            #facet.by = "Types",
            legend.title="MRscore (mean)",
            risk.table.col = "strata",
            surv.median.line = "hv",
            risk.table.y.text.col = T,
            risk.table.y.text = T)
 fit<-survfit(Surv(DSStime/365,DSS) ~ MRscore_group,ctype = 1,
              conf.type="log-log",
              data = data2)
 surv_median(fit)
 ggsurvplot(fit = fit,facet.by = "Stage_ct",
            data = data2,
            pval.method = T,
            risk.table = TRUE,
            pval = TRUE,
            ggtheme = theme_survminer(),
            palette = c("#ff9a3c","#283c63"),
            #facet.by = "Types",
            legend.title="MRscore (mean)",
            risk.table.col = "strata",
            surv.median.line = "hv",
            risk.table.y.text.col = T,
            risk.table.y.text = T)
 ##stage specific survival
 data3=subset(data2,Stage_ct=="Stage I")
 data3=subset(data2,Stage_ct=="Stage II")
 data3=subset(data2,Stage_ct=="Stage III")
 data3=subset(data2,Stage_ct=="Stage IV")
 data3=subset(data2,Types=="BLCA")
 fit<-survfit(Surv(OStime/365,OS) ~ MRscore_group,ctype = 1,
              conf.type="log-log",
              data = data3)
 surv_median(fit)
 ggsurvplot(fit = fit,facet.by = "Stage_ct",
            data = data3,
            pval.method = T,
            risk.table = TRUE,
            pval = TRUE,
            ggtheme = theme_survminer(),
            palette = c("red","blue"),
            #facet.by = "Types",
            legend.title="MRscore(mean)",
            risk.table.col = "strata",
            surv.median.line = "hv",
            risk.table.y.text.col = T,
            risk.table.y.text = T)
 ##glmnet continus variables cox
 library(glmnet)
 library(survival)
 d <- subset(meta,select = c(OStime,OS,MRscore,Age,Stage_ct))
 d<-d%>%subset(.,Age!="NA")%>%subset(.,Stage_ct!="NA")%>%subset(.,OStime!=0)
 x <- model.matrix( ~ MRscore + Age + Stage_ct, d)
 y <- Surv(d$OStime, d$OS)
 fit <- glmnet(x, y, family="cox", alpha=1)
 plot(fit, label=T)#顶端的横坐标应该是当前Lambda下非零变量的个数：
 #glmnet()返回的是一系列不同Lambda对应的值（一组模型），
 coef(fit)
 #需要user来选择一个Lambda，交叉验证是最常用挑选Lambda的方法
 #LASSO回归复杂度调整的程度由参数λ来控制，λ越大对变量较多的线性模型的惩罚力度就越大，从而最终获得一个变量较少的模型
 cv.fit <- cv.glmnet(x, y, family="cox", alpha=1,nfolds = 10)
 plot(cv.fit)
 #lambda.min（cross-validation误差最小）、
 #lambda.1se（cross-validation误差最小一个标准差内，模型最简单），对应图中两根竖线的地方
 #It includes the cross-validation curve (red dotted line), and upper and lower standard deviation（标准差） curves along
 #the lambda sequence (error bars). Two selected lambda’s are indicated by the vertical dotted lines (see below).

 #lambda.min is the value of Lambda that gives minimum mean cross-validated error. 
 coefficient <- coef(cv.fit, s = cv.fit$lambda.min)
 Active.index<-which(as.numeric(coefficient)!=0)
 Active.coefficients<-as.numeric(coefficient)[Active.index]
 multi_cox_res<-rownames(coefficient)[Active.coefficients]
 
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
library(forestmodel)
p<-forest_model(coxfit,factor_separate_line=F, 
                 format_options = list(colour= "black", 
                                       shape = 12, 
                                       text_size = 4, 
                                       banded = T), 
                 theme = theme_forest())
 p
 
 #cox indifferent stage
 cox_early_OS<- coxph(Surv(OStime, OS)~MRscore, 
             data = subset(meta,Stage=="Low_stage"))
 cox_late_OS<- coxph(Surv(OStime, OS)~MRscore, 
                   data = subset(data,Stage=="High_stage"))
 cox_early_PFI<- coxph(Surv(PFItime, PFI)~MRscore, 
                      data = subset(meta,Stage=="Low_stage"))
 cox_late_PFI<- coxph(Surv(PFItime, PFI)~MRscore, 
                     data = subset(meta,Stage=="High_stage"))
 cox_early_DSS<- coxph(Surv(DSStime, DSS)~MRscore, 
                       data = subset(meta,Stage=="Low_stage"))
 cox_late_DSS<- coxph(Surv(DSStime, DSS)~MRscore, 
                      data = subset(meta,Stage=="High_stage"))
 
 coef_early_OS<-summary(cox_early_OS)$coefficients%>%
   as.data.frame()%>%mutate(Group="early_OS")%>%
   mutate(Variables=rownames(summary(cox_early_OS)$coefficients))%>%
   cbind(.,summary(cox_early_OS)$Concordance)
 
 coef_late_OS<-summary(cox_late_OS)$coefficients%>%
   as.data.frame()%>%mutate(Group="late_OS")%>%
   mutate(Variables=rownames(summary(cox_late_OS)$coefficients))
 
 coef_early_PFI<-summary(cox_early_PFI)$coefficients%>%
   as.data.frame()%>%mutate(Group="early_PFI")%>%
   mutate(Variables=rownames(summary(cox_early_PFI)$coefficients))
 
 coef_late_PFI<-summary(cox_late_PFI)$coefficients%>%
   as.data.frame()%>%mutate(Group="late_PFI")%>%
   mutate(Variables=rownames(summary(cox_late_PFI)$coefficients))
 
 coef_early_DSS<-summary(cox_early_DSS)$coefficients%>%
   as.data.frame()%>%mutate(Group="early_DSS")%>%
   mutate(Variables=rownames(summary(cox_early_DSS)$coefficients))
 
 coef_late_DSS<-summary(cox_late_DSS)$coefficients%>%
   as.data.frame()%>%mutate(Group="late_DSS")%>%
   mutate(Variables=rownames(summary(cox_late_DSS)$coefficients))
 coef=bind_rows(list(coef_early_OS,coef_late_OS,
                     coef_early_PFI, coef_late_PFI,
                     coef_early_DSS, coef_late_DSS))
 data = subset(meta,Stage=="Low_stage")
 riskScore=predict(cox_early_OS,type="risk",newdata=data)
 risk=as.vector(ifelse(riskScore>1,"high","low"))
 COX_risk_result<-cbind(id=rownames(cbind(data[,8:9],riskScore,risk)),
                             cbind(data[,8:9],riskScore,risk))
 write.csv(COX_risk_result,"../dataset/TCGA_cox/MRscore_COX_riskinEarlystage.csv")
 rt=read.csv("../dataset/TCGA_cox/MRscore_COX_riskinEarlystage.csv",check.names=F,row.names=1)
 roc=survivalROC(Stime=rt$OStime, status=rt$OS, marker = rt$riskScore, 
                 predict.time =3, method="KM")
 pdf(file=paste0("../dataset/TCGA_cox/MRscore_COX_riskROC_early.pdf"))
 par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
 plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
      xlab="False positive rate", ylab="True positive rate",
      main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
      lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
 abline(0,1)
 dev.off()
 fit<- survfit(Surv(OStime/30, OS) ~ risk, data = rt)
 diff=survdiff(Surv(OStime/30, OS) ~risk,data = rt)
 pValue=1-pchisq(diff$chisq,df=1)
 pValue=round(pValue,14)
 plot(fit, lty = 2:3,col=c("red","blue"),xlab="time (months)",
      ylab="survival rate",
      main=paste("survival curve (p=", pValue ,")",sep=""),mark.time=T)
 legend("topright", c("Low risk", "High risk"), lty = 2:3, col=c("blue","red"))
 
for(i in c(1:6)){
    res_table[[i]]<-as.data.frame(cbind(cox_res[[i]]$coefficients,
                                   cox_res[[i]]$conf.int[,-c(1,2)]))
    }
 res_table<-as.data.frame(cbind(cox_res$coefficients,cox_res$conf.int[,-c(1,2)]))
 rownames(res_table)<-gsub("low","",rownames(mutil_res$coefficients))
 res_table<-subset(res_table,`Pr(>|z|)`<=0.05)
 
 
 
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
