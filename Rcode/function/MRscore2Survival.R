#' @ grid makeContent
#' 
#' 

  require(survival)
  require(survminer)
  require(survivalROC)
  require(forestplot)
  MRscore2Survival<- function(MRscore,meta,survivalTypes){
    survdata<-merge(MRscore,meta,by="sampleID")
    if(survivalTypes=="OS"){
      surv_cut <- surv_cutpoint(
        survdata,
        time = "OS",
        event = "Status",
        variables = c("MRscore")
      )
      summary(surv_cut)
      survdata$Group<-sample(c("High","Low"),nrow(MRscore),replace = T)
      survdata$Group<-ifelse(survdata$MRscore>=as.numeric(summary(surv_cut)[1]),"High","Low")
      fit<-survfit(Surv(OS,Status) ~ Group,
                   data = survdata)
      p<-ggsurvplot( fit,
                     data=survdata,
                     risk.table = TRUE,
                     pval = TRUE,
                     title = "OS",
                     palette = c("blue","red"),
                     #facet.by = "Efficacy",
                     legend.title="MIRscore",
                     risk.table.col = "strata",
                     surv.median.line = "hv",
                     risk.table.y.text.col = T,
                     risk.table.y.text = FALSE )
      p}
    else{
      surv_cut <- surv_cutpoint(
        survdata,
        time = "PFS",
        event = "Status.1",
        variables = c("MRscore")
      )
      summary(surv_cut)
      survdata$Group<-sample(c("High","Low"),nrow(MRscore),replace = T)
      survdata$Group<-ifelse(survdata$MRscore>=as.numeric(summary(surv_cut)[1]),"High","Low")
      fit<-survfit(Surv(PFS,Status.1) ~ Group,
                   data = survdata)
      p<-ggsurvplot( fit,
                     data=survdata,
                     risk.table = TRUE,
                     pval = TRUE,
                     title = "PFS",
                     palette = c("blue","red"),
                     #facet.by = "Efficacy",
                     legend.title="MIRscore",
                     risk.table.col = "strata",
                     surv.median.line = "hv",
                     risk.table.y.text.col = T,
                     risk.table.y.text = FALSE )
      
    }
    return(p,survdata)
  }
  MRscore2Sva(MRscore = MRscore,meta = meta,survivalTypes = "PFS")
  MRscore2Sva(MRscore = MRscore,meta = meta,survivalTypes = "OS")
  my_comparisons<-list(c("PD", "SD"), c("PD", "PR"), c("PR", "SD"))
  
  surv_cut <- surv_cutpoint(
    survdata,
    time = "PFS",
    event = "Status.1",
    variables = c("MRscore")
  )
  summary(surv_cut)
  plot(surv_cut)
  surv_cut <- surv_cutpoint(
    survdata,
    time = "OS",
    event = "Status",
    variables = c("MRscore")
  )
  summary(surv_cut)
  plot(surv_cut)
  survdata<-merge(MRscore,trans_meta,by="sampleID")
  survdata$MRscore_level<-sample(c("High","Low"),nrow(MRscore),replace = T)
  survdata$MRscore_level<-ifelse(survdata$MRscore>=as.numeric(summary(surv_cut)[1]),"High","Low")
  
  fun_to_plot <- function(data, group, variable) {
    p <- ggboxplot(data, x=group, y=variable,fill = group, 
                   palette = c("#00AFBB", "#E7B800", "#FC4E07"), 
                   add = "jitter", shape=group)+
      stat_compare_means(comparisons = my_comparisons)+
      stat_compare_means(label.y = 30)
    return(p)
  }