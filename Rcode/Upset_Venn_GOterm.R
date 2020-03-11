library(UpSetR)
library(tidyverse)
df<-read.csv("../dataset/GOterm_define/genesAnnotation.csv",header = T)
Fungus<-df[,-c(3,4,5,6,7)] %>% na.omit(.) %>% .[!duplicated(.$Symbol),]
antimic<-df[,-c(2,4,5,6,7)]%>% na.omit(.)%>%.[!duplicated(.$Symbol),]
Virus<-df[,-c(3,2,5,6,7)]%>% na.omit(.)%>%.[!duplicated(.$Symbol),]
Bac<-df[,-c(3,4,2,6,7)]%>% na.omit(.)%>%.[!duplicated(.$Symbol),]
InnateImmune<-df[,-c(3,4,5,2,7)]%>% na.omit(.)%>%.[!duplicated(.$Symbol),]
InflammatoryResponse<-df[,-c(3,4,5,6,2)]%>% na.omit(.)%>%.[!duplicated(.$Symbol),]
data<-list("Response to fungus" = Fungus$Symbol,
           "Response to bacterium"=Bac$Symbol,
           "Response to virus"=Virus$Symbol,
           "Inflammatory response"=InflammatoryResponse$Symbol,
           "Innate immune"=InnateImmune$Symbol,
           "Antimicrobial humoral response"=antimic$Symbol
           )
library(RColorBrewer)
colors()
p1<-upset(fromList(data),nsets = 6, 
          order.by = c("freq","degree"),matrix.color = "grey4",
          scale.intersections = "identity",
          sets.bar.color = c("purple","yellow2","royalblue3","peru","sienna","skyblue2"),point.size =4,
          set_size.numbers_size = F,
          set_size.scale_max = 950,
          set_size.show = T,
          mb.ratio = c(0.6, 0.4),
          text.scale = c(2, 1.5, 2, 1.5, 2, 1.5),
          sets.x.label = "GeneNumber",
          queries = list(list(query=intersects, params=list("Response to bacterium", "Inflammatory response"), color="purple", active=T), 
                                          list(query=intersects, params=list("Response to bacterium", "Innate immune"), color="purple", active=T), 
                                          list(query=intersects, params=list("Response to bacterium", "Inflammatory response","Innate immune"), color="purple", active=T),
                        list(query=intersects, params=list("Response to bacterium", "Innate immune","Antimicrobial humoral response"), color="purple", active=T),
                        list(query=intersects, params=list("Response to bacterium"), color="purple", active=T)
                        ))
p1




