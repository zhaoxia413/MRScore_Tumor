# ref:https://blog.csdn.net/weixin_43249695/article/details/102913239
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
cli<-fread("../dataset/TCGA_data/TCGA_clinical.csv")
colnames(cli)[1]<-"Tumor_Sample_ID"
rownames(cli)<-cli$Tumor_Sample_ID
MRS_expr<-fread("../dataset/TCGA_data/expr_MicroSignature_TMP_T.csv")
MRS_expr$Tumor_Sample_ID<-str_sub(MRS_expr$Tumor_Sample_ID,1,12)
head(cli$Tumor_Sample_ID)
MRS_expr<-MRS_expr[!duplicated(MRS_expr$Tumor_Sample_ID)]
rownames(MRS_expr)<-MRS_expr$Tumor_Sample_ID
head(MRS_expr$Tumor_Sample_ID)
cli_MRS_expr<-merge(cli[,1:4],MRS_expr,by="Tumor_Sample_ID")
cli_MRS_expr_cancers<-split.data.frame(cli_MRS_expr,cli_MRS_expr$type)
length(cli_MRS_expr_cancers)
cli_MRS_expr_cancers<-cli_MRS_expr_cancers[-which(sapply(cli_MRS_expr_cancers, nrow)<=100)]
sapply(cli_MRS_expr_cancers, nrow)
expr<-list()
surv<-list()
group_data<-list()
survival_dat<-list()
survival_dat_merge<-list()
checkGroup<-list()
covariates<-list()
for (i in seq_len(length(cli_MRS_expr_cancers))) {
  print(i)
  expr[[i]]<-cli_MRS_expr_cancers[[i]][,c(5:831)]
  rownames(expr[[i]])<-cli_MRS_expr_cancers[[i]]$Tumor_Sample_ID
  which(is.na(expr[[i]]))
  expr[[i]]<-as.data.frame(lapply(expr[[i]],as.numeric))
  surv[[i]]<-cli_MRS_expr_cancers[[i]][,c(1:4)]
  rownames(surv[[i]])<-cli_MRS_expr_cancers[[i]]$Tumor_Sample_ID
  group_data[[i]] <- apply(expr[[i]], 2 , function(gene){
    name <- colnames(gene)
    gene <- unlist(gene)
    group <- ifelse(gene >= median(gene), 'high', 'low')
    names(group) <- name
    return(group)
  })
  group_data[[i]] <- as.data.frame(group_data[[i]], stringsAsFactors = F)
  
  survival_dat[[i]] <- data.frame(row.names = rownames(surv[[i]]),status = surv[[i]]$OS,
                             time = surv[[i]]$OS.time,
                             stringsAsFactors = F)
  survival_dat_merge[[i]] <- cbind(survival_dat[[i]],group_data[[i]])
  names(survival_dat_merge)[i]<-levels(factor(cli$type))[i]
  checkGroup[[i]]<-apply(survival_dat_merge[[i]],2 , function(gene){
    facter_lenth<-length(levels(factor(gene)))
    check<-ifelse(facter_lenth==1,facter_lenth,"Yes")
    names(check)<-colnames(gene)
    return(check)
  })
  checkGroup[[i]] <- as.data.frame(checkGroup[[i]], stringsAsFactors = F)
  covariates[[i]] <- as.character(rownames(subset(checkGroup[[i]],checkGroup[[i]]!=1)))
}
univ_formulas<-list()
for (i in seq_len(length(covariates))) {
  univ_formulas[[i]] <- sapply(covariates[[i]],
                          function(x){
                            ##print(x)
                            as.formula(paste('Surv(time, status)~', x))
                          })
}

univ_models<-list()
univ_results<-list()
res_single<-list()
single_pick<-list()
univ_models<-list()
#21,22 error
#error test:
# for (i in 3:765) {
#  univ_models[[i]]<-coxph(univ_formulas[[23]][[i]], data = survival_dat_merge[[23]])
# }
sapply(univ_formulas, length)
for (i in seq_len(length(univ_formulas))) {
  print(paste0("i=",i))
  univ_models[[i]]<-list()
  for (j in 3:length(univ_formulas[[i]])){
    print(paste0("j=",j))
    univ_models[[i]][[j]]<-coxph(univ_formulas[[i]][[j]], data = survival_dat_merge[[i]])  
  }
  names(univ_models[[i]])<-univ_formulas[[i]][3:length(univ_formulas[[i]])]
}

for (i in seq_len(length(univ_models))) {
  print(i)
  univ_models[[i]]<-univ_models[[i]][-which(sapply(univ_models[[i]], is.null))]
  # Extract data
  univ_results[[i]] <- lapply(univ_models[[i]],
                              function(x){
                                x <- summary(x)
                                p.value <- signif(x$wald["pvalue"], digits = 2)
                                beta <- signif(x$coef[1], digits = 2)
                                HR <- signif(x$coef[2], digits = 2)
                                HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                                HR.confint.upper <- signif(x$conf.int[,"upper .95"], 2)
                                HR <- paste0(HR, " (",
                                             HR.confint.lower, "-", HR.confint.upper, ")")
                                res <- c(beta, HR, p.value)
                                names(res) <- c("coef", "HR (95% CI for HR)", "p.value")
                                return(res)
                              })
}
for (i in seq_len(length(univ_results))) {
  print(length(univ_results[[i]]))
}
for (i in seq_len(length(univ_results))) {
  print(i)
  res_single[[i]] <- as.data.frame(t(do.call(cbind, univ_results[[i]])))
  table(res_single[[i]]$p.value <= 0.05)
  res_single[[i]] <- res_single[[i]][res_single[[i]]$p.value <= 0.05, ]
  res_single[[i]] <- res_single[[i]][order(res_single[[i]]$p.value), ]
  rownames(res_single[[i]]) <- gsub('Surv.time..status....',"",rownames(res_single[[i]]))
  single_pick[[i]] <-rownames(res_single[[i]])
  names(single_pick)[i]<-levels(factor(cli$type))[i]
  write.csv(res_single[[i]],paste0("../dataset/COX_results/KM_res_single_TCGA_",levels(factor(cli$type))[i],".csv"))
}
# 3. Cox regressison
## multi
which(sapply(single_pick, is.null))
names(single_pick)
length(single_pick)
fmla<-list()
cox<-list()
mutil_res<-list()
mutil_res_table<-list()
multi_pick<-list()
final_res<-list()
riskScore<-list()
risk<-list()
multiCOX_risk_result<-list()
rt<-list()
roc<-list()
diff<-list()
pValue<-list()
fit<-list()
survival_dat_merge<-survival_dat_merge[names(single_pick)]
names(survival_dat_merge)
names(single_pick)
sapply(single_pick, length)
survival_dat_merge<-survival_dat_merge[which(sapply(single_pick, length)>10)]
single_pick<-single_pick[which(sapply(single_pick, length)>10)]
sapply(single_pick, length)
seq_len(length(single_pick))
for (i in seq_len(length(single_pick))) {
  print(i)
fmla[[i]] <- as.formula(paste0("Surv(time, status) ~",paste0(single_pick[[i]],collapse = '+')))
names(fmla)[i]<-names(single_pick)[i]
}
library(survcomp)
C_index<-list()

for (i in seq(fmla)) {
  tryCatch(
    {cox[[i]] <- coxph(fmla[[i]], data = survival_dat_merge[[i]])#构建多因素生存模型

  if(length(cox[[i]])==0){
    next
  }
mutil_res[[i]]<-summary(cox[[i]])
mutil_res_table[[i]]<-as.data.frame(cbind(mutil_res[[i]]$coefficients,mutil_res[[i]]$conf.int[,-c(1,2)]))
rownames(mutil_res_table[[i]])<-gsub("low","",rownames(mutil_res[[i]]$coefficients))
mutil_res_table[[i]]<-subset(mutil_res_table[[i]],`Pr(>|z|)`<=0.05)
multi_pick[[i]]<-rownames(mutil_res_table[[i]])
if(length(multi_pick[[i]])==0){
  next
}
print(paste0(names(fmla)[i],"past"))
write.csv(mutil_res_table[[i]],paste0("../dataset/COX_results/mutil_res_table_",names(fmla)[i],".csv"))
names(multi_pick)[i]<-names(single_pick)[i]
fmla[[i]] <- as.formula(paste0("Surv(time, status) ~",paste0(multi_pick[[i]],collapse = '+')))
cox[[i]] <- coxph(fmla[[i]], data = survival_dat_merge[[i]])#构建多因素生存模型
final_res[[i]]<-summary(cox[[i]])
cox[[i]]=step(cox[[i]],direction = "both")
riskScore[[i]]=predict(cox[[i]],type="risk",newdata=survival_dat_merge[[i]])
risk[[i]]=as.vector(ifelse(riskScore[[i]]>1,"high","low"))
names(risk)[i]<-names(multi_pick)[i]
multiCOX_risk_result[[i]]<-cbind(id=rownames(cbind(survival_dat_merge[[i]][,1:2],riskScore[[i]],risk[[i]])),
                                 cbind(survival_dat_merge[[i]][,1:2],riskScore[[i]],risk[[i]]))
write.csv(multiCOX_risk_result[[i]],paste0("../dataset/COX_results/mutation_model-multiCOX_risk_",names(survival_dat_merge)[i],".csv"))
rt[[i]]=read.csv(paste0("../dataset/COX_results/mutation_model-multiCOX_risk_",names(survival_dat_merge)[i],".csv"),check.names=F,row.names=1)
roc[[i]]=survivalROC(Stime=rt[[i]]$time, status=rt[[i]]$status, marker = rt[[i]]$riskScore, 
                     predict.time =3, method="KM")
pdf(file=paste0("../dataset/COX_results/mutation_model-multiCOX_riskROC_",names(survival_dat_merge)[i],".pdf"))
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
plot(roc[[i]]$FP, roc[[i]]$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",round(roc[[i]]$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()
# multiCOX_Survival Curve
diff[[i]]=survdiff(Surv(time, status) ~risk[[i]],data = rt[[i]])
pValue[[i]]=1-pchisq(diff[[i]]$chisq,df=1)
pValue[[i]]=round(pValue[[i]],14)
fit[[i]] <- survfit(Surv(time/30, status) ~ risk[[i]], data = rt[[i]])
names(cox)[i]=names(fmla)[i]
C_index[[i]]<-data.frame(concordance=summary(cox[[i]])$concordance)%>%mutate(Types=rep(names(cox)[i]),nrow(.))
pdf(file=paste0("../dataset/COX_results/mutation_model-multiCOX_survival_curve_",names(survival_dat_merge)[i],".pdf"))
plot(fit[[i]], lty = 2:3,col=c("red","blue"),xlab="time (months)",
      ylab="survival rate",
     main=paste("survival curve (p=", pValue[[i]] ,")",sep=""),mark.time=T)
legend("topright", c("Low risk", "High risk"), lty = 2:3, col=c("blue","red"))
dev.off()
    }, 
error=function(e) NA)}
Cindex=bind_rows(C_index)
Cindex$var=rep(c("C_index","se"),nrow(Cindex)/2)
write.csv(Cindex,"../dataset/COX_results/Cindex_cancerTypes.csv",row.names = F)

#test
sapply(cox, length)
#[1] 23 23 22 22  0 23  0 22 23  0  0 23  0 23 23 22
sapply(cox, length)

sapply(multi_pick, length)
#ACC BLCA BRCA CESC      COAD      ESCA  GBM           KIRC      LAML 
#8    8  122    0    0   31    0  111  136    0    0   21    0    5 
#LGG LIHC 
#1    0 
for (i in c(1:5)) {
 
  if (i==4){
   next
  }
  print(i)
}
#C-Index
length(cox)
for (i in seq(cox)) {
  summary(cox[[i]])
}
summary(cox[[1]])

library(survcomp)
library(survival)

fmla <- as.formula(paste0("Surv(time, status) ~",paste0(row.names(data.filtered)[1:10],collapse = '+')))
cox <- coxph(fmla, data = as.data.frame(t(data.filtered)))#构建多因素生存模型

cindex <- concordance.index(predict(cox),surv.time = time, surv.event = status ,method = "noether") cindex$c.index; cindex$lower; cindex$upper
