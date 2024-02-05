rm(list = ls())
setwd('~/eabenchmark/') 

library(tidyverse)

networks = c('funcoup','string')
ora = c('ease','fisher')
anubix = paste('anubix',networks,sep = '-')
neat = paste('neat',networks,sep = '-')
binox = paste('binox',networks,sep = '-')
netpea = paste('netpea',networks,sep = '-')
camera = c('camera','camera*')
methods = c(anubix,binox,neat,netpea,camera,ora,'roast','fgsea','gsa','gsea','padog','gsva','cepaORA','spia')

top_cutoff = c(10,20,40,'ALL')

message('READING POSITIVE-BENCHMARK')
positive_benchmark = read_tsv(file = 'results/stats/TPbenchmark.tsv.gz') %>%
  filter(Method %in% methods) %>% droplevels() %>% 
  replace(is.na(.),1) %>% 
  select(GEO,pathway_id,Pvalue,Method,label,TOP) %>%
  dplyr::rename(Pvalue=Pvalue,
                pathway_id = pathway_id) %>%
  add_column(benchmark = 'TP')

message('READING NEGATIVE-BENCHMARK')
class_balance = TRUE
if (class_balance){
  negative_benchmark = read_tsv(file = 'results/stats/FPbenchmark_balanced.tsv.gz') %>% 
    filter(Method %in% methods) %>% droplevels() %>% 
    replace(is.na(.),1) %>%
    select(GEO,pathway_id,Pvalue,Method,label,TOP) %>%
    dplyr::rename(Pvalue=Pvalue,
                  pathway_id = pathway_id) %>%
    add_column(benchmark = 'FP')
} else {
  negative_benchmark = read_tsv(file = 'results/stats/FPbenchmark.tsv.gz') %>% 
    replace(is.na(.),1) %>%
    filter(Method %in% methods) %>% droplevels() %>% 
    group_by(GEO,pathway_id,Method,Pvalue)%>%
    summarise(TOP = list(top_cutoff)) %>%
    unnest(cols = TOP) %>%
    ungroup() %>% 
    add_column(label = 1,
               benchmark = 'FP')
} 

################################################
## STATS OVERALL
benchmark_df = rbind(positive_benchmark,negative_benchmark)

stats_df = benchmark_df %>%
  group_by(Method,TOP) %>%
  ## Pvalue
  summarise(TP=as.double(sum(Pvalue<0.05 & label==1 & benchmark=='TP')),
            FN=as.double(sum(Pvalue>=0.05 & label==1 & benchmark=='TP')),
            FP=as.double(sum(Pvalue<0.05 & label==1 & benchmark=='FP')),
            TN=as.double(sum(Pvalue>=0.05 & label==1 & benchmark=='FP'))) %>%
  ungroup() %>%
  mutate(TPR = TP/(TP+FN),
         FPR = FP/(TN+FP) ) %>%
  mutate(GM = sqrt(TPR*(1-FPR))) 

if (class_balance){
  stats_df %>%
    write_tsv('results/stats/stats_balanced.tsv.gz')
} else {
  stats_df %>%
    write_tsv('results/stats.tsv.gz')
}