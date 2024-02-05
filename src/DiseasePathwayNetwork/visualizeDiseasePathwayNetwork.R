rm(list = ls())
setwd('~/eabenchmark/') 

library(tidyverse)
library(igraph)

##################################
### Target Disease Pathway Network
##################################
dataset_df = read_tsv(file='data/TPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv', 
                      col_types = 'fcfffiddiciccccc') %>% 
  select(GEO,TargetPathway) %>%
  droplevels()

top_cutoff = c(20)
pathway_network_df = lapply(top_cutoff, function(top){
  tmp_df = read_tsv(file = paste0('results/DiseasePathwayNetwork/DiseasePathwayNetwork-FDR<5e-2_top',top,'GOaverageBMA.tsv.gz'), col_types = cols())
  tmp_df$TOP = top
  return(tmp_df)
}) %>% data.table::rbindlist() %>%
  filter(TOP == '20') %>%
  rename(P1=A_name, P2=B_name) %>% 
  select(P1,P2,GO_BMA,TOP) %>% 
  droplevels()

# Full network
g = graph_from_data_frame(pathway_network_df,directed = F)
subg = igraph::simplify(g,remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = 'max')
vcount(subg); ecount(subg)

# Target network
pathway_network_df = pathway_network_df %>% 
  filter(P1 %in% dataset_df$TargetPathway)
g = graph_from_data_frame(pathway_network_df,directed = F)
## All pathways
subg = igraph::simplify(g,remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = 'max')
vcount(subg); ecount(subg)

kegg_df = read_tsv(file = 'data/KEGG/keggrest_hsa.tsv', col_types = 'ccccccc') %>%
  select(pathway_name,pathway_id,subclass) %>% 
  unique() %>%
  droplevels()

color_df = data.frame(
  subclass = c(unique(kegg_df %>%
                        filter(pathway_name %in% dataset_df$TargetPathway) %>% 
                        pull(subclass))),
  color = c(
    '#ffdd55',
    '#d90429',
    '#f7aef8',
    '#49b6ff',
    '#04a777',
    '#ff9e00'
  ))

kegg_df = kegg_df %>%
  left_join(color_df, by='subclass') %>% 
  replace(is.na(.), '#eaebed')
kegg_df


label_df = data.frame(TargetPathway = unique(dataset_df$TargetPathway)) %>% 
  rowwise() %>%
  mutate(subclass = kegg_df[kegg_df$pathway_name==TargetPathway,'subclass']) %>% 
  arrange(subclass) %>%
  ungroup() %>% 
  mutate(label = 1:nrow(.))

# E(subg)$width = scale01(E(subg)$GO_BMA)
E(subg)$width = E(subg)$GO_BMA
V(subg)$color = unlist(sapply(V(subg)$name,function(p) as.character(kegg_df[kegg_df$pathway_name==p,'color'])), use.names = FALSE)
V(subg)$frame.width = 2
V(subg)$size = unlist(sapply(V(subg)$name,function(p) ifelse(p %in% dataset_df$TargetPathway, 4, 2)), use.names = FALSE)
V(subg)$frame.color = rgb(red=0, green=0, blue=0, maxColorValue=255, alpha = 255)
V(subg)$frame.width = unlist(sapply(V(subg)$name,function(p) ifelse(p %in% dataset_df$TargetPathway, 2, 1.5)), use.names = FALSE)
# V(subg)$width = 1
V(subg)$label = unlist(sapply(V(subg)$name,function(p) ifelse(p %in% dataset_df$TargetPathway, as.character(label_df[label_df$TargetPathway==p,][['label']]), NA)), use.names = FALSE)
V(subg)$shape = 'square'


data.frame(p = V(subg)$name, color = V(subg)$color) %>%
  left_join(color_df, by='color') 

set.seed(17) 
#The LGL algorithm is meant for large, connected graphs. 
#Here you can also specify a root: a node that will be placed in the middle of the layout.

plot(subg,
     # vertex.label = NA,
     vertex.label.dist = 0,
     vertex.label.dist = -1,
     vertex.label.position = 'left',
     vertex.label.color = 'black',
     vertex.label.family = "sans",
     vertex.label.cex = 0.5, 
     vertex.label.font = 2,
     vertex.shape = "square",
     edge.label.cex = .5)

tmp_color_df = color_df %>% 
  filter(subclass %in% (kegg_df %>%
                          filter(pathway_name %in% V(subg)$name) %>%
                          pull(subclass))) %>% 
  arrange(as.character(subclass)) %>% droplevels()

tmp_color_df = rbind(tmp_color_df,data.frame(subclass='Not disease pathway',color='#eaebed'))

legend(-0.4,-0.6,bty = "n",
       cex = 0.85,
       legend=tmp_color_df$subclass,
       fill=tmp_color_df$color, border=T)
