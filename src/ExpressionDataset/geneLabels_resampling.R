rm(list = ls())
setwd('~/eabenchmark/') 

library(tidyverse)
library(SummarizedExperiment)

n_sampling = 30
indir = 'data/TPbenchmark/'
outdir = 'data/FPbenchmark/'

## FUNCOUP mapping file
uniprot_mapping_df = 
  read_tsv('data/mapping/UP000005640_9606.idmapping.gz', 
           col_names = c('Uniprot_id','db','db_id'), 
           col_types = 'cfc') %>%
  mutate(db_id=gsub('hsa:','',db_id)) 

uniprot_mapping_df = uniprot_mapping_df %>%
  filter(db=='GeneID') %>%
  select(Uniprot_id,db_id) %>% 
  dplyr::rename(GeneID=db_id) %>%
  inner_join(uniprot_mapping_df %>%
               filter(db=='Ensembl'), by = 'Uniprot_id', multiple = "all") %>% 
  select(GeneID,db_id) %>% 
  dplyr::rename(Ensembl=db_id) %>%
  unique()

## REMOVE ENSEMBL --> GENEID for >1 times
uniprot_mapping_df = 
  uniprot_mapping_df %>% 
  filter(Ensembl %in% 
           (uniprot_mapping_df %>% 
              group_by(Ensembl) %>%
              summarise(TOT=n()) %>% 
              filter(TOT==1) %>%
              pull(Ensembl)))

string_mapping_df = 
  read_tsv('data/mapping/9606.protein.aliases.v11.5.txt.gz', 
           col_names = c('string_protein_id','alias','source'), 
           col_types = 'ccf', comment = '#') %>%
  filter(grepl('Gene_ID|GeneID',source, perl = T)) %>%
  mutate(string_protein_id = gsub('9606\\.','',string_protein_id)) %>%
  dplyr::rename(Ensembl_PRO=string_protein_id,
                GeneID=alias) %>%
  select(GeneID,Ensembl_PRO) %>% 
  unique()

## REMOVE ENSEMBL_PRO --> GENEID for >1 times
string_mapping_df = 
  string_mapping_df %>% 
  filter(Ensembl_PRO %in% 
           (string_mapping_df %>% 
              group_by(Ensembl_PRO) %>%
              summarise(TOT=n()) %>% 
              filter(TOT==1) %>%
              pull(Ensembl_PRO)))

## Join Uniprot and STRING mapping where no Ensembl/Ensembl_PRO id maps to multiple GeneID
## The opposite can still happen, meaning GeneID can map to multiple Ensembl 
mapping_df = 
  uniprot_mapping_df %>% 
  inner_join(string_mapping_df, by = 'GeneID', multiple = 'all')

mapping_df %>% 
  select(GeneID,Ensembl) %>% 
  unique() %>% 
  group_by(Ensembl) %>%
  summarise(TOT=n()) %>% 
  arrange(-TOT) %>% head

mapping_df %>% 
  select(GeneID,Ensembl_PRO) %>% 
  unique() %>% 
  group_by(Ensembl_PRO) %>%
  summarise(TOT=n()) %>% 
  arrange(-TOT) %>% head

human_genes = unique(mapping_df$GeneID)

#################################################
## EXPRESSION DATA
dataset_df = read_tsv(file='data/TPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv', 
                      col_types = 'fcfffiddiciccccci')

count = 1; tot = nrow(dataset_df)
out = apply(dataset_df,1,function(x){
  message(paste0(x[['GEO']],' ', count,'/',tot)); count <<- count +1
  
  se = readRDS(paste0(indir,x[['DataType']],'/',x[['GEO']],'.Rds'))
  dataset_genes = rownames(se)
  
  set.seed(x[['seed.number']])
  seeds_perm = sample(10000,n_sampling,replace = F)
  if(length(dataset_genes)>length(human_genes)){
    permuted_se = 
      as.data.frame(
        sapply(
          seq(n_sampling), function(i){
            set.seed(seeds_perm[i])
            return(sample(dataset_genes,nrow(se),replace = F))
          })
      )
    
  } else {
    permuted_se = 
      as.data.frame(
        sapply(
          seq(n_sampling), function(i){
            set.seed(seeds_perm[i])
            return(sample(human_genes,nrow(se),replace = F))
          })
      )
  }
  
  colnames(permuted_se) = paste0('P',seq(n_sampling))
  write_tsv(permuted_se, file = paste0(outdir,x[['GEO']],'.tsv.gz'))
})

### RUN extract_deg.R --> extract DEG according to n_sampling just computed