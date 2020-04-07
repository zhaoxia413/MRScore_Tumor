library(tidyverse)
library(igraph)
library(ggraph)
library(reshape2)
library(data.table)

data <- read.csv("test_MRgeneset.csv",header = T)
# most popular programming languages from TIOBE Index (Nov. 2019) found in data
# (only languages with position <= 16 are considered)
popular_languages <- c("C1","C2","C3")
colnames(data)[1]
# number of packages to display
number_of_pkgs <- 300
# find largest packages written in popular languages
top_packages <- data %>%
  filter(language %in% popular_languages) %>%
  group_by(pkg_name) %>%
  summarize(total_code = sum(code)) %>%
  arrange(desc(total_code)) %>%
  head(number_of_pkgs) %>%
  select(pkg_name, total_code)

# all popular languages per package
top_languages_per_pkg <- data %>%
  arrange(pkg_name, desc(code)) %>%
  group_by(pkg_name) %>%
  mutate(
    main = row_number() == 1, # main language of package should be opaque
    total_code = sum(code)
  ) %>%
  ungroup() %>%
  select(language, pkg_name, code, total_code, main)

# only following languages found in given packages
(top_languages <- top_languages_per_pkg %>%
    pull(language) %>%
    unique %>%
    sort)

top_language_colors <- c(
  '#efb306',
  '#eb990c',
  '#000000'
)

names(top_language_colors) <- c(
  "C1","C2","C3"
)

edges1 <- top_languages_per_pkg %>%
  transmute(from = language, to = pkg_name, total_code = code, main)

edges2 <- top_languages_per_pkg %>%
  count(language, wt = code, name = 'total_code') %>%
  transmute(
    from = '',
    to = language,
    total_code,
    main = TRUE
  )

edges <- bind_rows(edges1, edges2)

vertices1 <- top_languages_per_pkg %>%
  filter(main) %>%
  transmute(
    node = pkg_name, language, total_code, level = 1
  )

vertices2 <- edges2 %>%
  transmute(
    node = to, language = to, total_code, level = 2
  )

vertices3 <- tibble(
  node = '', language = NA, total_code = 0, level = 3
)

vertices <- bind_rows(vertices1, vertices2, vertices3) %>%
  mutate(
    radius = total_code**(1.8), # scaling circles
    language = factor(language, names(top_language_colors))
  ) %>%
  arrange(level, language, node)
graph <- graph_from_data_frame(edges, vertices = vertices)

# create custom layout by updating existing circle layout
layout <- create_layout(graph, layout = 'circle')

outer_circle <- layout %>%
  filter(level == 1) %>%
  mutate(language = factor(language, names(top_language_colors))) %>%
  arrange(language, desc(name)) %>%
  mutate(
    x = cos((row_number() - 1) / number_of_terms * 2 * pi),
    y = sin((row_number() - 1) / number_of_terms * 2 * pi)
  )

# positioning circle centers manually by specifying polar coords
angles <- c(3, 43, 50,0)
radii <- c(0.8, 0.38, 0.5,0)
centers <- tibble(
  x = radii * cos(angles / 180 * pi),
  y = radii * sin(angles / 180 * pi)
)
inner_circle <- bind_cols(centers, select(filter(layout, level != 1), -x, -y))

layout[] <- bind_rows(outer_circle, inner_circle) %>%
  arrange(ggraph.index)

ggraph(layout) +
  geom_edge_diagonal(
    aes(edge_color = node1.language, edge_alpha = as.factor(main)),
    edge_width = 0.3, show.legend = FALSE
  ) +
  geom_node_point(
    aes(size = radius, color = language),
    alpha = 0.6, show.legend = FALSE
  ) +
  geom_node_text(
    aes(
      x = 1.0175 * x,
      y = 1.0175 * y,
      label = name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = !(name %in% top_languages)
    ),
    size = 2, hjust = 'outward', family = 'Oswald'
  ) +
  geom_node_text(
    aes(
      x = x,
      y = y,
      label = name,
      filter = name %in% top_languages
    ),
    size = 6, hjust = 0.5, family = 'Oswald'
  ) +
  geom_node_text(
    aes(
      x = x,
      y = y - 0.045,
      label = ifelse(
        total_code > 1000,
        format(total_code, big.mark = ','),
        total_code
      ),
      filter = name %in% top_languages
    ),
    size = 3, hjust = 0.5, family = 'Oswald'
  ) +
  scale_edge_color_manual(values = top_language_colors) +
  scale_color_manual(values = top_language_colors) +
  scale_size_area(max_size = 150) +
  scale_edge_alpha_manual(values = c(0.15, 1)) +
  coord_fixed() +
  labs(
    title = 'LOC of Popular Programming Languages in 300 CRAN Packages',
    subtitle = 'considered are largest CRAN packages written in one (or more) of top 16 programming languages from TIOBE Index (Nov. 2019)',
    caption = '#tidytuesday 46|2019 spren9er'
  ) +
  theme_void() +
  theme(
    text = element_text(family = 'Oswald'),
    legend.position = c(0.645, 0.51),
    plot.title = element_text(
      face = 'bold', hjust = 0.5, size = 20, margin = margin(t = 45, b = 3)
    ),
    plot.subtitle = element_text(
      face = 'plain', hjust = 0.5, size = 13, margin = margin(t = 5, b = 3)),
    plot.caption = element_text(
      face = 'plain', color = '#dedede', size = 8, hjust = 1,
      margin = margin(b = 20)
    )
  )

ggsave(
  'images/tidytuesday_201946_cran_packages.png',
  width = 12, height = 12.5, dpi = 300
)
library(data.table)
library(tidyverse)
options(stringsAsFactors = F)
df<-read.csv("../dataset/MRgenesets/MRgeneSet.csv",header = T)
head(df)
length(levels(factor(df$Annotated.Term)))#465
length(levels(factor(df$Symbol)))#2519
signature<-fread("../dataset/TCGA_data/annotationRow1.csv")%>%as.data.frame()
head(signature)
df1<-filter(df,df$Symbol%in%signature$Gene)
ConsensusGenes<-c("CXCL10","CXCL9","CXCL1","CXCL3","CYBA","PYCARD",
                  "GBP5","OAS2","IL1B","STAT1","OAS3","DEFA1","DEFA3","CXCL2",
                  "IL6","TLR4","GBP2","IL18","CD14","LY96","IFNAR2","CD86","CASP1"
                  ,"TLR2","GBP3","IRF5")
ConsensusGenesAnno<-df1[which(df1$Symbol%in%ConsensusGenes),]
write.csv(ConsensusGenesAnno,"ConsensusGenesAnno.csv",row.names = F)
ConsensusGenesAnno<-ConsensusGenesAnno[,c(2,5)]%>%group_by(Symbol,Annotated.Term)%>%summarise(n())
ConsensusGenesAnno<-ConsensusGenesAnno[,c(1,2)]
ConsensusGenesAnno<-ConsensusGenesAnno%>%group_by(Symbol,Annotated.Term)%>%summarise(n())
set1<-c("CXCL10","CXCL9","CXCL1","CXCL3","CYBA","PYCARD",
        "GBP5","OAS2","IL1B","STAT1","OAS3")
set2<-c("DEFA1","DEFA3","CXCL2",
        "IL6","TLR4")
set3<-c("GBP2","IL18","CD14","LY96","IFNAR2","CD86","CASP1"
        ,"TLR2","GBP3","IRF5")
Con1<-ConsensusGenesAnno[which(ConsensusGenesAnno$Symbol%in%set1),]
Con2<-ConsensusGenesAnno[which(ConsensusGenesAnno$Symbol%in%set2),]
Con3<-ConsensusGenesAnno[which(ConsensusGenesAnno$Symbol%in%set3),]
plot1<-Con1%>%group_by(Annotated.Term)%>%summarise(Count=n())
plot2<-Con2%>%group_by(Annotated.Term)%>%summarise(Count=n())
plot3<-Con3%>%group_by(Annotated.Term)%>%summarise(Count=n())
topgene1<-top_n(plot1,5,wt = plot1$Count)%>%mutate(Class=rep("Clster1",nrow(.)))
topgene2<-top_n(plot2,5,wt = plot2$Count)%>%mutate(Class=rep("Clster2",nrow(.)))
topgene3<-top_n(plot3,5,wt = plot3$Count)%>%mutate(Class=rep("Clster3",nrow(.)))
data<-bind_rows(topgene1,topgene2,topgene3)
library(ggsci)
library(ggthemes)
ggplot(data,aes(Class,Count,fill=Annotated.Term))+
  geom_col()+
  theme_few()+
  scale_fill_d3()

length(levels(factor(df1$Symbol)))#1202
length(levels(factor(df1$Annotated.Term)))
length(levels(factor(df1$Evidence)))#13
head(df1)
df1%>%group_by(Evidence)%>%summarise(n())
# Evidence `n()`
# <fct>    <int>
# 1 IBA       2361 Inferred from biological aspect of ancestor
# 2 IC          22
# 3 IDA       1019 Inferred from direct assay
# 4 IEA       1412 Inferred from electronic annotation
# 5 IEP         62
# 6 IGI        243
# 7 IMP       1654 Inferred from mutant phenotype
# 8 IPI         13
# 9 ISO       2042 Inferred from sequence orthology
# 10 ISS         12
# 11 NAS         20
# 12 RCA          2
# 13 TAS         72
colnames(df1)
df2<-df1%>%group_by(Goterm,Symbol,Chr,Annotated.Term)%>%summarise(Dupli=n())
length(levels(factor(df2$Symbol)))
df3<-df2%>%group_by(Goterm,Annotated.Term)%>%summarise(Dupli=n())
datalist<-split.data.frame(df3,f=df3$Goterm,drop = F)
mytop<-function(x){
  top_n(x,5,wt = x$Dupli)
}
datalist1<-lapply(datalist, mytop)
dftop10<-bind_rows(datalist1)
head(dftop10)
library(ggalluvial)
library(reshape2)
library(tidyverse)
head(dftop10)
geneAnnotation<-df1[,c(2,5)]
geneAnnotation<-geneAnnotation%>%group_by(Symbol,Annotated.Term)%>%summarise(n())%>%.[,-3]
data<-merge(dftop10,geneAnnotation,by="Annotated.Term")
colnames(data)[4]<-"Gene"
MRDEGs<-fread("../dataset/TCGA_data/microSignature_DEGtag.csv")[,c(1,4,6)]%>%as.data.frame()
data1<-merge(data,MRDEGs,by="Gene",all = T)
levels(factor(data1$Goterm))
data1<-subset(data1,Annotated.Term!="NA")
data1$Goterm<-gsub("DefenseResponse","Def_R",data1$Goterm)
data1$Goterm<-gsub("ImmuneResponse","Im_R",data1$Goterm)
data1$Goterm<-gsub("Response2ExBioticStimulus","R_ExBio",data1$Goterm)
library(RColorBrewer)
col31<-c("#303841","#D72323","#377F5B","#375B7F","#F2FCFC","#f0027f",
         "#FAF8DE","#666666","#BDF1F6","#023782","#5e4fa2","#F1C40F",
         "#ff7f00","#cab2d6","#240041","#ffff99","#0E3BF0","#a65628",
         "#f781bf","#808FA6","#2EB872","#F0FFE1","#F33535","#011F4E",
         "#82B269","#D3C13E","#3F9DCD","#014E1F","#AFFFDF","#3D002E",
         "#3A554A")
ggplot(data = data1,
       aes(axis1 = Goterm, axis2 = Annotated.Term,axis3=Types,
           weight =  Dupli)) +
  scale_x_discrete(limits = c("Dupli", "Annotated.Term","Types"), expand = c(.01, .05)) +
  geom_alluvium(aes(fill = Types)) +
  geom_stratum() + geom_text(stat = "stratum", label.strata = T) +
  #theme_minimal() +
  #theme_tropical()
  #theme_matrix()
  theme_bw(base_size = 8)+
  scale_fill_manual(values = col31)+
  guides(fill=NULL)

suppressPackageStartupMessages( require(easyalluvial))
library(ggthemes)
termnamme<-paste0("Term",seq(1,19,1))
termnamme[1:9]<-c("Term01","Term02","Term03","Term04",
                  "Term05","Term06","Term07","Term08","Term09")
termAnno<-data.frame(Annotated.Term=levels(factor(data1$Annotated.Term)),Annotated.Term_T=termnamme)
write.csv(termAnno,"termAnno.csv",row.names = F)
data2<-merge(data1,termAnno,by="Annotated.Term")
data2=data1
data2<-filter(data1,Annotated.Term!="immune response")
head(data2)
data2$Regulation<-data2$Log2FC
data2$Regulation<-ifelse(data2$Regulation>0,"Up",ifelse(data2$Regulation<0,"Down","A1"))
alluvial_wide(data2[,c(3,2,5,7,6)],
              fill_by = 'first_variable',
               stratum_labels = T,
               #col_vector_value =col31,
               col_vector_flow = col31[c(16,17,19)],
               colorful_fill_variable_stratum = T,
               stratum_label_size = 4,
               stratum_width = 1/4
               , order_levels = c('8','6','4')
               )+
  theme_few(base_size = 12)
par(mfrow=c(10,1));
par(mar=c(0.1,0.1,2,0.1), xaxs="i", yaxs="i")
barplot(rep(1,times=10),col=col31[c(1:3,5,6,13:17)],border=cm.colors(10),
        axes=FALSE, main="cm.colors"); box()
library(pheatmap)
mat<-data.frame(row.names = data2$Gene,log2FC=data2$Log2FC)

# find largest packages written in popular languages
number_of_terms <- 300
df2$Dupli<-rep(1,nrow(df2))
df2<-df2[df2$Annotated.Term%in%dftop20$Annotated.Term,]
length(levels(factor(df2$Annotated.Term)))
top_terms <- df2 %>%
  group_by(Annotated.Term) %>%
  summarize(total_gene = n()) %>%
  arrange(desc(total_gene)) %>%
  head(number_of_terms) %>%
  select(Annotated.Term, total_gene)

# all popular languages per package
top_terms_per_GO <- df2 %>%
  arrange(Annotated.Term, desc(Dupli)) %>%
  group_by(Annotated.Term) %>%
  mutate(
    main = row_number() == 1, # main language of package should be opaque
    total_gene = sum(Dupli)
  ) %>%
  ungroup() %>%
  select(Goterm, Annotated.Term, Dupli, total_gene, main)

# only following languages found in given packages
(top_terms <- top_terms_per_GO %>%
    pull(Goterm) %>%
    unique %>%
    sort)

top_GO_colors <- c(
  '#efb306',
  '#eb990c',
  '#000000'
)

names(top_GO_colors) <- c(
  "DefenseResponse","ImmuneResponse","Response2ExBioticStimulus"
)

edges1 <- top_terms_per_GO %>%
  transmute(from = Goterm, to = Annotated.Term, total_gene = Dupli, main)

edges2 <- top_terms_per_GO %>%
  count(Goterm, wt = Dupli, name = 'total_gene') %>%
  transmute(
    from = '',
    to = Goterm,
    total_gene,
    main = TRUE
  )

edges <- bind_rows(edges1, edges2)

vertices1 <- top_terms_per_GO %>%
  filter(main) %>%
  transmute(
    node = Annotated.Term, Goterm, total_gene, level = 1
  )

vertices2 <- edges2 %>%
  transmute(
    node = to, Goterm = to, total_gene, level = 2
  )

vertices3 <- tibble(
  node = '', Goterm = NA, total_gene = 0, level = 3
)

vertices <- bind_rows(vertices1, vertices2, vertices3) %>%
  mutate(
    radius = total_gene**(1.8), # scaling circles
    Goterm = factor(Goterm, names(top_GO_colors))
  ) %>%
  arrange(level, Goterm, node)
graph <- graph_from_data_frame(edges, vertices = vertices)

# create custom layout by updating existing circle layout
layout <- create_layout(graph, layout = 'circle')

outer_circle <- layout %>%
  filter(level == 1) %>%
  mutate(Goterm = factor(Goterm, names(top_GO_colors))) %>%
  arrange(Goterm, desc(name)) %>%
  mutate(
    x = cos((row_number() - 1) / number_of_terms * 2 * pi),
    y = sin((row_number() - 1) / number_of_terms * 2 * pi)
  )

# positioning circle centers manually by specifying polar coords
angles <- c(3, 43, 50,0)
radii <- c(0.8, 0.38, 0.5,0)
centers <- tibble(
  x = radii * cos(angles / 180 * pi),
  y = radii * sin(angles / 180 * pi)
)
inner_circle <- bind_cols(centers, select(filter(layout, level != 1), -x, -y))

layout[] <- bind_rows(outer_circle, inner_circle) %>%
  arrange(ggraph.index)

ggraph(layout) +
  geom_edge_diagonal(
    aes(edge_alpha = as.factor(main)),
    edge_width = 0.3, show.legend = FALSE
  ) +
  geom_node_point(
    aes(size = radius, color = Goterm),
    alpha = 0.6, show.legend = FALSE
  ) +
  geom_node_text(
    aes(
      x = 1.0175 * x,
      y = 1.0175 * y,
      label = name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = !(name %in% top_terms)
    ),
    size = 2, hjust = 'outward', family = 'Oswald'
  ) +
  geom_node_text(
    aes(
      x = x,
      y = y,
      label = name,
      filter = name %in% top_terms
    ),
    size = 6, hjust = 0.5, family = 'Oswald'
  ) +
  geom_node_text(
    aes(
      x = x,
      y = y - 0.045,
      label = ifelse(
        total_gene > 1000,
        format(total_gene, big.mark = ','),
        total_gene
      ),
      filter = name %in% top_terms
    ),
    size = 3, hjust = 0.5, family = 'Oswald'
  ) +
  scale_edge_color_manual(values = top_GO_colors) +
  scale_color_manual(values = top_GO_colors) +
  scale_size_area(max_size = 150) +
  scale_edge_alpha_manual(values = c(0.15, 1)) +
  coord_fixed() +
  labs(
    title = 'GO term',
    subtitle = 'GOterm',
    caption = '#goterm'
  ) +
  theme_void() +
  theme(
    text = element_text(family = 'Oswald'),
    legend.position = c(0.645, 0.51),
    plot.title = element_text(
      face = 'bold', hjust = 0.5, size = 20, margin = margin(t = 45, b = 3)
    ),
    plot.subtitle = element_text(
      face = 'plain', hjust = 0.5, size = 13, margin = margin(t = 5, b = 3)),
    plot.caption = element_text(
      face = 'plain', color = '#dedede', size = 8, hjust = 1,
      margin = margin(b = 20)
    )
  )

ggsave(
  'images/tidytuesday_201946_cran_packages.png',
  width = 12, height = 12.5, dpi = 300
)

