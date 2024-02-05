rm(list = ls())
setwd('~/eabenchmark/') 

networks = c('funcoup','string')
anubix = paste('anubix',networks,sep = '-')
neat = paste('neat',networks,sep = '-')
binox = paste('binox',networks,sep = '-')
netpea = paste('netpea',networks,sep = '-')
net_methods = c(anubix,binox,neat,netpea)

library(tidyverse)
source('src/others/loadNetwork.R')
library(gridExtra)
library(grid)
library(scales)
library(ggpubr)
library(ggridges)
library(parallel)
library(gridExtra)
library(grid)
range01 = function(x){(x-min(x))/(max(x)-min(x))}

# Plot style
convert_name=function(m){
  if(m == 'anubix-funcoup') return('ANUBIX-FunCoup')
  if(m == 'anubix-string') return('ANUBIX-STRING')
  if(m == 'binox-funcoup') return('BinoX-FunCoup')
  if(m == 'binox-string') return('BinoX-STRING')
  if(m == 'neat-funcoup') return('NEAT-FunCoup')
  if(m == 'neat-string') return('NEAT-STRING')
  if(m == 'netpea-funcoup') return('netPEA-FunCoup')
  if(m == 'netpea-string') return('netPEA-STRING')
}

plain_gg_theme=theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.border = element_rect(linetype = "solid", colour = "black", size=0),
        axis.line.y.right = element_line(colour = 'white', size = 0),
        axis.line.x.top = element_line(colour = 'white', size = 0),
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        plot.title=element_text(size=13, face='bold'),
        axis.title=element_text(size=13, face='bold'),
        strip.text=element_text(size=13, face='bold'),
        strip.background = element_rect(colour="white", fill="white"),
        legend.justification=c('center', 'top'),
        legend.box.just='right',
        legend.margin=margin(6, 6, 6, 6),
        legend.title= element_text(size=10, face='bold'),
        legend.text= element_text(size=10))

dataset_df = 
  read_tsv(file='data/TPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv', 
           col_types = 'fcfffiddiciccccc') %>%
  arrange(TargetPathway)

kegg_df = read_tsv(file = 'data/KEGG/keggrest_hsa.tsv', col_types = 'ffcccff')
kegg_df = kegg_df %>%
  filter(pathway_id %in%
           as.character(
             kegg_df %>%
               select(pathway_id,entrezgene_id) %>%
               unique() %>%
               group_by(pathway_id) %>%
               summarise(n=n()) %>%
               filter(n>14) %>%
               droplevels() %>%
               pull(pathway_id))) %>%
  droplevels() %>%
  dplyr::rename(GeneID=entrezgene_id)

kegg_list = tapply(kegg_df$GeneID, kegg_df$pathway_id, function(x) as.character(x))
kegg_list = lapply(kegg_list, function(x) unique(x))
geneset_names = names(kegg_list)

kegg_df = kegg_df %>% 
  select(pathway_id,pathway_name,class,subclass) %>% 
  unique() %>% 
  droplevels()

### Pathway Network properties
networks = c('funcoup','string')
network_colors = c('FunCoup'='#cc2936','STRING' ='#08415c')

net_property_df = lapply(networks, function(network){
  graph = load_network(network) 
  d = degree(graph)
  avg_degree = mean(d)
  message(paste0(ifelse(network=='funcoup','FunCoup','STRING')), ' avg. degree = ',round(avg_degree,2))
  kegg_df = read_tsv(paste0('data/KEGG/', network,'_kegg.tsv')) %>%
    dplyr::select(Pathway,Ensembl) %>%
    dplyr::rename(pathway_id = Pathway) %>%
    dplyr::filter(Ensembl %in% V(graph)$name) %>%
    unique()

  min_d = data.frame(d) %>%  
    arrange(-d) %>%
    head((nrow(.)/100)*20) %>% 
    tail(1) %>% pull(d)
  
  top20p_df = data.frame(d) %>%  
    arrange(-d) %>% 
    filter(d>=min_d) %>%
    mutate(V = rownames(.))
  rownames(top20p_df) = NULL
  
  tmp_df = lapply(unique(kegg_df$pathway_id),function(pathway){
    g = kegg_df %>%
      dplyr::filter(pathway_id==pathway) %>%
      pull(Ensembl)
    sub_graph = igraph::induced_subgraph(graph, which(V(graph)$name %in% g))
    sub_d = degree(sub_graph)
    data.frame(pathway_id = pathway,
               size_pathway = vcount(sub_graph),
               num_top20Hubs_net = nrow(top20p_df),
               num_top20Hubs = sum(sapply(g,function(v){
                 v %in% top20p_df$V
               })),
               num_components = igraph::count_components(sub_graph),
               mean_clustering = igraph::transitivity(sub_graph, type = 'average'),
               mean_degree = mean(d[g],na.rm = T),
               median_degree = median(d[g], na.rm = T),
               max_degree = max(d[g], na.rm = T),
               num_intralinks = sum(sub_d[g]),
               ratio_intralinks = sum(sub_d[g]) / sum(d[g]),
               ratio_full = (sum(sub_d[g])/2) / (vcount(sub_graph)**2)
    )
  }) %>% data.table::rbindlist()
  tmp_df$Network = network
  return(tmp_df)
}) %>% data.table::rbindlist()

net_property_df$Network = ifelse(net_property_df$Network == 'funcoup','FunCoup','STRING')
net_property_df$density_top20Hubs = net_property_df$num_top20Hubs / net_property_df$size_pathway #DENSITY of hubs
net_property_df$ratio_top20Hubs = net_property_df$num_top20Hubs / net_property_df$num_top20Hubs_net #FRACTION of hubs
net_property_df %>% 
  select(Network,num_top20Hubs_net) %>% 
  unique()


negative_benchmark = read_tsv('results/stats/FPbenchmark.tsv.gz')

negative_benchmark_df=negative_benchmark %>%
  filter(Method %in% net_methods &
           pathway_id %in% kegg_df$pathway_id
  ) %>%
  droplevels() %>%
  group_by(Method,pathway_id) %>%
  summarise(FPR=sum(Pvalue<0.05)/length(Pvalue)) %>%
  ungroup() %>%
  pivot_wider(names_from=Method, values_from=FPR) 

pathway_bias = negative_benchmark_df %>%
  tidyr::pivot_longer(!c(pathway_id), names_to = 'Method', values_to = 'FPR') %>%
  mutate(Network = ifelse(grepl('funcoup',Method),'FunCoup','STRING')) %>%
  left_join(net_property_df, by=c('pathway_id','Network'))

pathway_bias$Method = sapply(pathway_bias$Method, function(x) convert_name(x))
pathway_bias$Method = gsub('-[A-Z|a-z]*','',pathway_bias$Method)
pathway_bias$Method = as.factor(pathway_bias$Method)
pathway_bias$Network = as.factor(pathway_bias$Network)

## FPR vs Hubs
pathway_bias_tile = pathway_bias %>% 
  group_by(Method,Network) %>%
  replace(is.na(.), 0) %>% 
  summarise("Size" = list(as.numeric(unlist(cor.test(size_pathway,FPR, method = 'spearman')[c('estimate','p.value')]))),
            "Fraction of intralinks" = list(as.numeric(unlist(cor.test(ratio_intralinks,FPR, method = 'spearman')[c('estimate','p.value')]))),
            "Density of links" = list(as.numeric(unlist(cor.test(ratio_full,FPR, method = 'spearman')[c('estimate','p.value')]))),
            "Max degree" = list(as.numeric(unlist(cor.test(max_degree,FPR, method = 'spearman')[c('estimate','p.value')]))),
            "Median degree" = list(as.numeric(unlist(cor.test(median_degree,FPR, method = 'spearman')[c('estimate','p.value')]))),
            "Number of CC" = list(as.numeric(unlist(cor.test(num_components,FPR, method = 'spearman')[c('estimate','p.value')]))),
            "Mean clustering coefficient" = list(as.numeric(unlist(cor.test(mean_clustering,FPR, method = 'spearman')[c('estimate','p.value')]))),
            "Density of hubs" = list(as.numeric(unlist(cor.test(density_top20Hubs,FPR, method = 'spearman')[c('estimate','p.value')]))), 
            "Fraction of hubs" =  list(as.numeric(unlist(cor.test(ratio_top20Hubs,FPR, method = 'spearman')[c('estimate','p.value')]))),
            ) %>% 
  ungroup() %>% 
  pivot_longer(names_to = 'type', values_to = 'Correlation', cols = c("Size","Fraction of intralinks","Density of links","Max degree","Median degree","Number of CC","Mean clustering coefficient","Fraction of hubs","Density of hubs")) %>% 
  tidyr::unnest_wider(Correlation, names_sep = "_") %>% 
  dplyr::rename(Rho = Correlation_1, Pvalue = Correlation_2)

pathway_bias_tile %>% 
  write_tsv('notebooks/tables/RhoOfFPRvsNetProperties_NEA.tsv')

string_pathway_bias_tile = pathway_bias_tile %>% 
  filter(Network == 'STRING') %>%
  select(-c(Rho)) %>% 
  mutate(Pvalue=format(Pvalue,scientific = TRUE,digits = 2)) %>% 
  pivot_wider(names_from=Method, values_from=Pvalue)
funcoup_pathway_bias_tile = pathway_bias_tile %>% 
  filter(Network == 'FunCoup') %>%
  select(-c(Rho)) %>% 
  mutate(Pvalue=format(Pvalue,scientific = TRUE,digits = 2)) %>% 
  pivot_wider(names_from=Method, values_from=Pvalue)

rbind(funcoup_pathway_bias_tile,string_pathway_bias_tile) %>% 
  write_tsv('notebooks/tables/PvalueOfFPRvsNetProperties_NEA.tsv')

pathway_bias_tile = read_tsv('notebooks/tables/RhoOfFPRvsNetProperties_NEA.tsv')

pathway_bias_tile$type = factor(pathway_bias_tile$type, levels = c("Density of hubs", "Density of links","Fraction of hubs","Fraction of intralinks","Max degree","Median degree","Mean clustering coefficient","Number of CC","Size"))
pathway_bias_tile$cor_label = paste0(round(pathway_bias_tile$Rho,2),ifelse(pathway_bias_tile$Pvalue>0.05,'', '*'))

ggplot(pathway_bias_tile, aes(x=type, y=Method, fill=Rho)) +
  geom_tile(color="white",
            lwd=1,
            linetype=1) +
  geom_text(aes(label = cor_label)) +
  xlab('') + ylab('') +
    scale_fill_gradientn(expression("r"['s']),
                         colors=c('#cc2936','white','#08415c'),
                       limits=c(-1,1),
                       breaks=c(-1,-0.75,-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)) + 
  facet_wrap(~Network, scales="free") +
  plain_gg_theme + 
  theme(axis.text.y=element_text(size=13, face='bold',colour="black"),
        axis.text.x=element_text(size=13, face='bold', angle=45, hjust=1,colour="black"),
        strip.text=element_text(size=13, face='bold',colour="black"),
        legend.justification=c("center"),
        legend.position="top", legend.key.width=unit(2, 'cm'),
        legend.title=element_text(size = 15, vjust=0.75,face = 'bold'),
        legend.text = element_text(size=9))

##########
net_property_df = net_property_df %>%
  replace(is.na(.), 0)

funcoup_pathway = net_property_df %>% 
  filter(Network=='FunCoup') %>% 
  pull(pathway_id)
funcoup_prop_df = data.frame(net_property_df %>% 
  filter(Network=='FunCoup') %>% 
  select(-c(pathway_id,Network,num_top20Hubs_net)) %>% 
  scale()) %>%
  mutate(pathway_id = funcoup_pathway,
         Network = 'FunCoup')

string_pathway = net_property_df %>% 
  filter(Network=='STRING') %>% 
  pull(pathway_id)
string_prop_df = data.frame(net_property_df %>% 
                               filter(Network=='STRING') %>% 
                               select(-c(pathway_id,Network,num_top20Hubs_net)) %>% 
                               scale()) %>%
  mutate(pathway_id = string_pathway,
         Network = 'STRING')

net_property_df = rbind(funcoup_prop_df,string_prop_df) %>% 
  rename("Size" = 'size_pathway',
         "Fraction of intralinks" = 'ratio_intralinks',
         "Density of links" = 'ratio_full',
         "Max degree" = 'max_degree',
         "Mean degree" = 'mean_degree',
         "Median degree" = 'median_degree',
         "Number of CC" = 'num_components',
         "Mean clustering coefficient" = 'mean_clustering',
         "Density of hubs" = 'density_top20Hubs',
         "Fraction of hubs" = 'ratio_top20Hubs'
  ) %>% 
  select(c("Network","Fraction of hubs","Fraction of intralinks","Density of hubs", "Density of links","Max degree","Median degree","Mean clustering coefficient","Number of CC","Size")) %>%
  pivot_longer(names_to = 'type', values_to = 'Value', cols = c("Size","Fraction of intralinks","Density of links","Max degree","Median degree","Number of CC","Mean clustering coefficient","Fraction of hubs","Density of hubs"))

ggplot(net_property_df, aes(x = Value, y = type, fill = Network)) +
  geom_density_ridges() +
  ylab('') +
  xlab('') +
  facet_wrap(~Network) + 
  scale_fill_manual(values=list('FunCoup'='#014f86','STRING'='#80b5db')) +
  scale_x_continuous("z-score", n.breaks = 10) +
  plain_gg_theme + 
  theme(axis.title.y=element_text(size=13,face='bold',colour="black", vjust=5),
        axis.text.y=element_text(size=12,face='bold',colour="black"),
        axis.text.x=element_text(size=13, face='bold',colour="black"),
        legend.position='none',
        plot.margin=unit(c(1,1,1,1), 'cm'))


net_methods = c(anubix,binox,neat,netpea)
df=negative_benchmark %>%
  filter(Method %in% net_methods &
           pathway_id %in% kegg_df$pathway_id
  ) %>%
  droplevels() %>%
  group_by(Method,pathway_id) %>%
  summarise(n=sum(Pvalue<0.05)/length(Pvalue)) %>%
  ungroup() %>%
  pivot_wider(names_from=Method, values_from=n) 

df2_funcoup = apply(df[,net_methods[grepl('funcoup',net_methods)]],2,function(x){
  return(rank(x, ties.method = 'average'))
}) %>% as.data.frame() %>%
  add_column(pathway_id = df$pathway_id, .before = 'anubix-funcoup')

df2_funcoup$mean_bias=apply(df2_funcoup,1,function(x){
  return(mean(as.numeric(x[-1],na.rm=T)))
}) 

df2_funcoup_top10 = df2_funcoup %>%
  left_join(kegg_df, by='pathway_id') %>% 
  arrange(-mean_bias) %>% head(10)

df2_funcoup_top10_pathwayname = unique(df2_funcoup_top10$pathway_name)

df2_funcoup_top10_tile = df2_funcoup_top10 %>% 
  left_join(funcoup_prop_df,by = 'pathway_id') %>% 
  rename("Size" = 'size_pathway',
         "Fraction of intralinks" = 'ratio_intralinks',
         "Density of links" = 'ratio_full',
         "Max degree" = 'max_degree',
         "Mean degree" = 'mean_degree',
         "Median degree" = 'median_degree',
         "Number of CC" = 'num_components',
         "Mean clustering coefficient" = 'mean_clustering',
         "Density of hubs" = 'density_top20Hubs',
         "Fraction of hubs" = 'ratio_top20Hubs'
  ) %>% 
  select(c("pathway_name", "Network","Fraction of hubs","Fraction of intralinks","Density of hubs", "Density of links","Max degree","Median degree","Mean clustering coefficient","Number of CC","Size")) %>%
  pivot_longer(names_to = 'type', values_to = 'Value', cols = c("Size","Fraction of intralinks","Density of links","Max degree","Median degree","Number of CC","Mean clustering coefficient","Fraction of hubs","Density of hubs"))

## STRING
df2_string = apply(df[,net_methods[grepl('string',net_methods)]],2,function(x){
  return(rank(x, ties.method = 'average'))
}) %>% as.data.frame() %>%
  add_column(pathway_id = df$pathway_id, .before = 'anubix-string')

df2_string$mean_bias=apply(df2_string,1,function(x){
  return(mean(as.numeric(x[-1],na.rm=T)))
}) 

df2_string_top10 = df2_string %>%
  left_join(kegg_df, by='pathway_id') %>% 
  arrange(-mean_bias) %>% head(10)

df2_string_top10_pathwayname = unique(df2_string_top10$pathway_name)

df2_string_top10_tile = df2_string_top10 %>% 
  left_join(string_prop_df,by = 'pathway_id') %>% 
  rename("Size" = 'size_pathway',
         "Fraction of intralinks" = 'ratio_intralinks',
         "Density of links" = 'ratio_full',
         "Max degree" = 'max_degree',
         "Mean degree" = 'mean_degree',
         "Median degree" = 'median_degree',
         "Number of CC" = 'num_components',
         "Mean clustering coefficient" = 'mean_clustering',
         "Density of hubs" = 'density_top20Hubs',
         "Fraction of hubs" = 'ratio_top20Hubs'
  ) %>% 
  select(c("pathway_name", "Network","Fraction of hubs","Fraction of intralinks","Density of hubs", "Density of links","Max degree","Median degree","Mean clustering coefficient","Number of CC","Size")) %>%
  pivot_longer(names_to = 'type', values_to = 'Value', cols = c("Size","Fraction of intralinks","Density of links","Max degree","Median degree","Number of CC","Mean clustering coefficient","Fraction of hubs","Density of hubs"))


tmp_df_top10 = rbind(df2_funcoup_top10_tile,df2_string_top10_tile)
tmp_df_top10$pathway_name = as.character(tmp_df_top10$pathway_name)
tmp_df_top10$pathway_name = ifelse(tmp_df_top10$pathway_name=='Chemical carcinogenesis - reactive oxygen species', 'Chemical carcinogenesis (ROS)',tmp_df_top10$pathway_name )
tmp_df_top10$pathway_name = ifelse(tmp_df_top10$pathway_name=='Pathways of neurodegeneration - multiple diseases', 'Pathways of neurodegeneration',tmp_df_top10$pathway_name )

top_biased_pathways = unique(c(df2_funcoup_top10_pathwayname,df2_string_top10_pathwayname))
intersect(df2_funcoup_top10_pathwayname,df2_string_top10_pathwayname)
tmp_df_top10$type = factor(tmp_df_top10$type, levels = c("Density of hubs", "Density of links","Fraction of hubs","Fraction of intralinks","Max degree","Median degree","Mean clustering coefficient","Number of CC","Size"))
tmp_df_top10$pathway_name = factor(tmp_df_top10$pathway_name, levels=top_biased_pathways)

gmt1 = ggplot(tmp_df_top10, aes(x=pathway_name, y=type, fill=Value)) +
  geom_tile(color="white",
            lwd=1,
            linetype=1) +
  xlab('') + ylab('') +
  scale_fill_gradientn('z-score',
    colors=c('#cc2936','white','#08415c'),
     values = scales::rescale(c(round(min(tmp_df_top10$Value)), 0, round(max(tmp_df_top10$Value)))),
     breaks=seq(-2,6,1)) + 
  facet_wrap(~Network, scales="free") +
  plain_gg_theme + 
  theme(axis.text.y=element_text(size=8, face='bold',colour="black"),
        axis.text.x=element_text(size=8, face='bold', angle=45, hjust=1,colour="black"),
        strip.text=element_text(size=8, face='bold',colour="black"),
        legend.justification=c("center"),        
        legend.position="top", legend.key.width=unit(1.5, 'cm'),
        legend.title=element_text(vjust=0.75,face = 'bold'),
        legend.text = element_text(size=9))

write_tsv(tmp_df_top10,file = 'notebooks/tables/top10fpr_networkPropertiesRange01.tsv')


## FPR vs pathways
df_funcoup_top10 = df[,c('pathway_id',net_methods[grepl('funcoup',net_methods)])]
df2_funcoup_top10_tile = df_funcoup_top10 %>% 
  filter(pathway_id %in% df2_funcoup_top10$pathway_id) %>% droplevels() %>% 
  left_join(kegg_df,by = 'pathway_id') %>% 
  pivot_longer(!c(pathway_id,pathway_name,class,subclass), names_to='Method', values_to='FPR') %>%
  mutate(Network = 'FunCoup')


df_string_top10 = df[,c('pathway_id',net_methods[grepl('string',net_methods)])]
df2_string_top10_tile = df_string_top10 %>% 
  filter(pathway_id %in% df2_string_top10$pathway_id) %>% droplevels() %>% 
  left_join(kegg_df,by = 'pathway_id') %>% 
  pivot_longer(!c(pathway_id,pathway_name,class,subclass), names_to='Method', values_to='FPR') %>%
  mutate(Network = 'STRING')

df_10 = rbind(df2_funcoup_top10_tile,df2_string_top10_tile)
df_10$Method=factor(df_10$Method, levels=sort(net_methods))
df_10$Method=sapply(df_10$Method, function(x) convert_name(x))
df_10$pathway_name = factor(df_10$pathway_name, levels=top_biased_pathways)

gmt2 = ggplot(df_10, aes(x=pathway_name, y=Method, fill=FPR)) +
  geom_tile(color="white",
            lwd=1,
            linetype=1) +
  geom_text(aes(label=round(FPR,2)), color="black", size=3) +
  xlab('') + ylab('') +
  scale_fill_gradientn(colors=c('white','#cc2936'),
                      values = scales::rescale(c(0.05, max(df_10$FPR))),
                       limits = c(0.05, max(df_10$FPR)),
                       breaks=seq(0.05,max(df_10$FPR),0.05)) + 
  facet_wrap(~Network, scales="free") +
  plain_gg_theme + 
  theme(axis.text.y=element_text(size=8, face='bold',colour="black"),
        axis.text.x=element_text(size=8, face='bold', angle=45, hjust=1,colour="black"),
        strip.text=element_text(size=8, face='bold',colour="black"),
        legend.title=element_text(vjust=0.75,face = 'bold'),
        legend.justification=c("center"),        
        legend.position="top", legend.key.width=unit(1.5, 'cm'),
        legend.text = element_text(size=9))

write_tsv(df_10,file = 'notebooks/tables/top10fpr_FunCoupANDString.tsv')

grid.arrange(gmt2,gmt1,ncol=1,heights=c(0.45,0.55))
# Add labels to the panels
grid.text("A", x = 0.05, y = 0.97, gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("B", x = 0.05, y = 0.53, gp = gpar(fontsize = 14, fontface = "bold"))