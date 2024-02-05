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

########################
## POSITIVE BENCHMARK
positive_benchmark = 
  read_tsv(file = 'results/stats/TPbenchmark.tsv.gz', col_types = 'fcdcff') %>% 
  replace(is.na(.),1)

top_cutoff = c(10,20,40,'ALL')
disease_network = lapply(top_cutoff, function(top){
  tmp_df = read_tsv(file = paste0('results/DiseasePathwayNetwork/DiseasePathwayNetwork-FDR<5e-2_top',top,'GOaverageBMA.tsv.gz'), col_types = cols())
  tmp_df$TOP = top
  return(tmp_df)
}) %>% data.table::rbindlist()

rank_df = lapply(methods, function(m){
  message(paste0('\n',m))
  top_cutoff = c(10,20,40,'ALL')
  m_positive_benchmark = positive_benchmark %>%
    filter(Method == m) %>% droplevels()
  
  tmp_dataset_df = dataset_df %>%
    dplyr::select(GEO,TargetPathway,pathway_id) %>% 
    droplevels()
  tmp_disease_network = disease_network
  
  compute_rank = 
    function(geo){
      geo = as.character(geo)
      print(geo)
      library(tidyverse)
      out_df = lapply(top_cutoff, function(t){
        target_p =
         tmp_dataset_df %>%
           filter(GEO==geo) %>%
           select(TargetPathway,pathway_id)
        
        ##target p
        set.seed(12345)
        tp_rank = m_positive_benchmark %>%
         filter(GEO == geo &
                  TOP == t ) %>%
         dplyr::slice(sample(1:n())) %>%
         ## Solve ties with average method
         mutate(rank = rank(Pvalue, ties.method = 'average')) %>%
         filter(pathway_id==target_p$pathway_id) %>%
         pull(rank)
        
        tmp_pathwayome_df = tmp_disease_network %>%
         filter(A == target_p[['pathway_id']] &
                  TOP == t) %>%
         add_row(A = target_p[['pathway_id']],
                 B = target_p[['pathway_id']],
                 A_name = target_p[['TargetPathway']],
                 B_name = target_p[['TargetPathway']],
                 TOP = t)
        
        ##target p + related pathways
        set.seed(12345)
        tmp_positive_benchmark = m_positive_benchmark %>%
         filter(GEO == geo &
                  TOP == t ) %>%
         dplyr::slice(sample(1:n())) %>%
         mutate(rank = rank(Pvalue, ties.method = 'average')) %>%
         left_join(tmp_pathwayome_df %>%
                     select(B) %>%
                     dplyr::rename(pathway_id=B), by = 'pathway_id')
        
        all_ranks = tmp_positive_benchmark %>%
         filter(label==1) %>%
         pull(rank)
        
        return(data.frame(GEO = geo,
                         TOP = t,
                         rank_target = tp_rank/318,
                         rank_all_median = median(all_ranks,na.rm=T)/318,
                         rank_all_mean = mean(all_ranks,na.rm=T)/318
                         )
       )
      })  %>% data.table::rbindlist()
      return(out_df)
    }

  pb = utils::txtProgressBar(min=0, max=length(unique(dataset_df$GEO)), style = 3)
  progress = function(n) utils::setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  cl = parallel::makeCluster(10)
  doSNOW::registerDoSNOW(cl)
  boot = foreach::foreach(i = unique(dataset_df$GEO), .options.snow = opts)
  out2_df = foreach::`%dopar%`(boot, compute_rank(i)) %>% 
    data.table::rbindlist() %>%
    add_column(Method = m)
  parallel::stopCluster(cl)
 
  return(out2_df)
}) %>% data.table::rbindlist()

write_tsv(rank_df, file = 'results/stats/Ranking(AVG).tsv.gz')
