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
                       aes(x=0.5,
                           0.17 * max(log(data[,y]+1), na.rm = TRUE),
                           label = text)) + 
      geom_text(data = resCor2, color="red",size=4,
                aes(x=0.5,
                    0.12 * max(log(data[,y]+1), na.rm = TRUE),
                    label = text))+
      theme_igray(base_size = 12)
  }
  return(p) 
}
