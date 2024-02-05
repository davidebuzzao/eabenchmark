rm(list = ls())
setwd('~/eabenchmark/') 

library(tidyverse)
library(parallel)

networks = c('funcoup','string')
ora = c('ease','fisher')
anubix = paste('anubix',networks,sep = '-')
neat = paste('neat',networks,sep = '-')
binox = paste('binox',networks,sep = '-')
netpea = paste('netpea',networks,sep = '-')
camera = c('camera','camera*')
methods = c(anubix,binox,neat,netpea,camera,ora,'roast','fgsea','gsa','gsea','padog','gsva','cepaORA','spia')


### Load data
kegg_df = read_tsv(file = 'data/KEGG/keggrest_hsa.tsv', col_types = 'cccccff')
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

dataset_df = 
  read_tsv(file='data/TPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv', 
           col_types = 'fcfffiddiciccccci') %>%
  left_join(kegg_df %>%
              dplyr::rename(TargetPathway='pathway_name'),by = 'TargetPathway')

dataset_df %>% 
  group_by(TargetPathway) %>% summarise(tot=n()) %>% arrange(-tot)

diseases = dataset_df %>%
  group_by(TargetPathway) %>% 
  summarise(tot=n()) %>%
  filter(tot>1) %>% 
  arrange(-tot) %>% 
  pull(TargetPathway)

disease_dataset_df = lapply(diseases, function(d){
  geos = as.character(dataset_df %>%
    filter(TargetPathway == d) %>% 
    pull(GEO))
  
  geos_df = lapply(utils::combn(unique(geos), 2, simplify = F),function(x){
    data.frame(A=x[1],B=x[2])
  } ) %>% data.table::rbindlist()
  geos_df$TargetPathway = d
  return(geos_df)
}) %>% data.table::rbindlist()

measure_reproducibility = function(m){
  library(tidyverse)
  tmp_benchmark = positive_benchmark %>% 
    filter(Method==m & Pvalue<0.05)
  
  tmp_reproducibility_df = apply(disease_dataset_df,1, function(x){
    A_pred = tmp_benchmark %>% 
      filter(GEO==x[['A']]) %>% 
      pull(pathway_id)
    B_pred = tmp_benchmark %>% 
      filter(GEO==x[['B']]) %>% 
      pull(pathway_id)
    
    data.frame(A = x[['A']],
               B = x[['B']],
               TargetPathway = x[['TargetPathway']],
               type = c('Jaccard index','Szymkiewicz–Simpson coefficient'),
               val = c(length(intersect(A_pred,B_pred))/length(union(A_pred,B_pred)),
                      length(intersect(A_pred,B_pred))/min(length(A_pred),length(B_pred))))
  }) %>% data.table::rbindlist()
  
  tmp_reproducibility_df$Method = m
  return(tmp_reproducibility_df)
}

########################
## POSITIVE BENCHMARK
positive_benchmark = 
  read_tsv(file = 'results/stats/TPbenchmark.tsv.gz', col_types = 'fcdcff') %>% 
  filter(TOP==20) %>%
  replace(is.na(.),1) %>% 
  droplevels()

methods = unique(positive_benchmark$Method)

pb = utils::txtProgressBar(min=0, max=length(methods), style = 3)
progress = function(n) utils::setTxtProgressBar(pb, n)
opts = list(progress = progress)
cl = parallel::makeCluster(8)
doSNOW::registerDoSNOW(cl)
boot = foreach::foreach(i = methods, .options.snow = opts)
reproducibility_df = foreach::`%dopar%`(boot, measure_reproducibility(i)) %>% 
  data.table::rbindlist()
parallel::stopCluster(cl)

write_tsv(reproducibility_df,'notebooks/tables/Reproducibility_JI&SS.tsv')

## PLOTTING
networks=c('funcoup','string')
ora = c('ease','fisher')
anubix = paste('anubix',networks,sep = '-')
neat = paste('neat',networks,sep = '-')
binox = paste('binox',networks,sep = '-')
netpea = paste('netpea',networks,sep = '-')
camera = c('camera','camera*')
methods = c(ora,camera,'roast','fgsea','gsa','gsea','padog','gsva','cepaORA','spia',anubix,binox,neat,netpea)


plain_gg_theme=theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.spacing = unit(2, "lines"),
        panel.border = element_rect(linetype = "solid", colour = "black", size=0),
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        plot.title=element_text(size=13, face='bold'),
        axis.title=element_text(size=13, face='bold'),
        strip.text=element_text(size=13, face='bold'),
        strip.background = element_rect(colour="white", fill="white"),
        # legend.box.just='right',
        legend.margin=margin(6, 6, 6, 6),
        legend.title= element_text(size=10, face='bold'),
        legend.text= element_text(size=10))

assign_category=function(m){
  if(m == 'fisher') return('Overlap-based')
  if(m == 'ease') return('Overlap-based')
  if(m == 'camera') return('Per-gene scoring (CP)')
  if(m == 'camera*') return('Per-gene scoring (CP)')
  if(m == 'gsea') return('Per-gene scoring (CP)')
  if(m == 'fgsea') return('Per-gene scoring (CP)')
  if(m == 'gsa') return('Per-gene scoring (CP)')
  if(m == 'padog') return('Per-gene scoring (CP)')
  if(m == 'gsva') return('Per-gene scoring (CP)')
  if(m == 'roast') return('Per-gene scoring (SC)')
  if(m == 'anubix-funcoup') return('FunCoup-based method')
  if(m == 'anubix-string') return('STRING-based method')
  if(m == 'binox-funcoup') return('FunCoup-based method')
  if(m == 'binox-string') return('STRING-based method')
  if(m == 'neat-funcoup') return('FunCoup-based method')
  if(m == 'neat-string') return('STRING-based method')
  if(m == 'netpea-funcoup') return('FunCoup-based method')
  if(m == 'netpea-string') return('STRING-based method')
  if(m == 'cepaORA') return('Pathway Topology-based')
  if(m == 'spia') return('Pathway Topology-based')
}

convert_name=function(m){
  if(m == 'fisher') return('Fisher')
  if(m == 'ease') return('EASE')
  if(m == 'camera') return('CAMERA_fixed')
  if(m == 'camera*') return('CAMERA_flex')
  if(m == 'gsea') return('GSEA')
  if(m == 'fgsea') return('fGSEA')
  if(m == 'gsa') return('GSA')
  if(m == 'padog') return('PADOG')
  if(m == 'roast') return('ROAST')
  if(m == 'gsva') return('GSVA')
  if(m == 'anubix-funcoup') return('ANUBIX-FunCoup')
  if(m == 'anubix-string') return('ANUBIX-STRING')
  if(m == 'binox-funcoup') return('BinoX-FunCoup')
  if(m == 'binox-string') return('BinoX-STRING')
  if(m == 'neat-funcoup') return('NEAT-FunCoup')
  if(m == 'neat-string') return('NEAT-STRING')
  if(m == 'netpea-funcoup') return('netPEA-FunCoup')
  if(m == 'netpea-string') return('netPEA-STRING')
  if(m == 'cepaORA') return('CePa-ORA')
  if(m == 'spia') return('SPIA')
}

reproducibility_df = read_tsv('notebooks/tables/Reproducibility_JI&SS.tsv')
reproducibility_avg_df = reproducibility_df %>% 
  group_by(Method,TargetPathway,type) %>%
  summarise(val = median(val,na.rm=T)) %>%
  ungroup() %>%
  replace(is.na(.),0)

disease_dataset_df %>% 
  group_by(TargetPathway) %>% 
  summarise(tot=n()) %>% arrange(-tot)

reproducibility_avg_df$TargetPathway = factor(reproducibility_avg_df$TargetPathway, levels = diseases)
order_m = as.character(sapply(methods, function(x) convert_name(x)))
reproducibility_avg_df$Method=sapply(reproducibility_avg_df$Method, function(x) convert_name(x))
reproducibility_avg_df$Method = factor(reproducibility_avg_df$Method, levels = order_m)

ggplot(reproducibility_avg_df, aes(x=Method, y=TargetPathway, fill=val)) +
  geom_tile(color="white",
            lwd=1,
            linetype=1) +
  xlab('') + ylab('') +
  scale_fill_gradientn(colors=c('white','#cc2936'),
                       limits=c(0,1),
                       breaks=seq(0,1,0.25)) + 
  facet_wrap(~type, scales="fixed") +
  plain_gg_theme + 
  theme(axis.text.y=element_text(size=13, face='bold',colour="black"),
        axis.text.x=element_text(size=13, face='bold', angle=45, hjust=1,colour="black"),
        strip.text=element_text(size=13, face='bold',colour="black"),
        legend.justification=c("center"),
        legend.position="top", legend.key.width=unit(2, 'cm'),
        legend.title=element_blank(),
        legend.text = element_text(size=9))

