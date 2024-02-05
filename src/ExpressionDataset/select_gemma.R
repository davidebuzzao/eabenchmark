rm(list = ls())
setwd('~/eabenchmark/') 

library(gemma.R)
library(tidyverse)
library(ggrepel)
library(SummarizedExperiment)

suppressMessages(library(parallel))

## GEO metadata were used to annotate our expression datasets.
## They can be fetched from GEO
ma.df = readRDS('data/geo/GEOmetadb_human_array.Rds')
rseq.df = readRDS('data/geo/GEOmetadb_human_rnaseq.Rds')

organism = 'hsa'
kegg_disease = read_tsv(paste0('data/KEGG/keggrest_',organism,'.tsv')) %>%
  filter(class=='Human Diseases') %>%
  pull(pathway_name) %>% 
  unique()

# Specify the directory path you want to create
dir_path = paste0("data/gemma/kegg_disease/",organism)

# Check if the directory exists, and create it if not
if (!file.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

fetch_info = function(pathway){
    message(pathway)
    if (file.exists(file = paste0('data/gemma/kegg_disease/',organism,'/',pathway,'.tsv'))){
      lookup = read_tsv(paste0('data/gemma/kegg_disease/',organism,'/',pathway,'.tsv'))
    } else {
      
      query = pathway
      ## Download data
      lookup =  tryCatch(
        {as_tibble(searchDatasets(query, limit = 100, taxon = ifelse(organism=='hsa','human','mouse'))) },
        error=function(cond) {
          return(NULL)})
      if(is.null(lookup)){return()}
      if (length(lookup)==0){return()}
      if (dim(lookup)[1]>0){
        
        lookup = lookup[ lookup$geeq.batchCorrected == TRUE |
                           (lookup$geeq.batchCorrected == FALSE & lookup$ee.batchEffect %in% c('No batch effect was detected')),] %>% 
          select(ee.ShortName, ee.Name, ee.ID, ee.Accession, ee.Samples, ee.batchEffect) 
        
        if (dim(lookup)[1]>0){
          lookup$query = query
          write_tsv(lookup, file = paste0('data/gemma/kegg_disease/',organism,'/',pathway,'.tsv'))
        }
      }
    }
    return(
      data.frame(
        pathway=pathway,
        datasets = ifelse(dim(lookup)[1]>0,dim(lookup)[1],0))
    )
  } 

ncores = 8
pb = utils::txtProgressBar(min=0, max=length(kegg_disease), style = 3)
progress = function(n) utils::setTxtProgressBar(pb, n)
opts = list(progress = progress)
cl = parallel::makeCluster(ncores)
doSNOW::registerDoSNOW(cl)
boot = foreach::foreach(i = kegg_disease, .options.snow = opts)
stats_df = foreach::`%dopar%`(boot, fetch_info(i)) %>% 
  data.table::rbindlist() %>% 
  unique()
parallel::stopCluster(cl)

stats_df = stats_df[stats_df$datasets>0,]
write_tsv(stats_df, paste0('data/gemma/kegg_disease/',organism,'_stats.tsv'))


#### Collect data
stats_df = read_tsv(paste0('data/gemma/kegg_disease/',organism,'_stats.tsv'),col_types = 'fi')
stats_df = stats_df %>% arrange(-datasets)

gemma_df = lapply( stats_df$pathway, function(pathway){
  read_tsv(paste0('data/gemma/kegg_disease/',organism,'/',pathway,'.tsv'),col_types = 'fccficf')
}) %>% data.table::rbindlist()

## Filter away not interesting ones
gemma_df = gemma_df[!grepl('treated|treatment|time',gemma_df$ee.Name, perl = T, ignore.case = T),]
gemma_df = gemma_df[grepl('GSE',gemma_df$ee.Accession),]
gemma_df$data.type = NA
gemma_df = gemma_df[gemma_df$ee.Samples>=8,]
gemma_df = unique(gemma_df)
write_tsv(gemma_df, paste0('data/gemma/kegg_disease/',organism,'/summary.tsv'))

## fetch gemma annotation
gemma_df = read_tsv( paste0('data/gemma/kegg_disease/',organism,'/summary.tsv'))
get_gemma_annotation = function(dataset){
  library(dplyr)
  library(gemma.R)
  getDatasetAnnotations(dataset) %>% 
    mutate(ee.shortName = dataset)
}
ncores = 8
pb = utils::txtProgressBar(min=0, max=length(unique(gemma_df$ee.ShortName)), style = 3)
progress = function(n) utils::setTxtProgressBar(pb, n)
opts = list(progress = progress)
cl = parallel::makeCluster(ncores)
doSNOW::registerDoSNOW(cl)
boot = foreach::foreach(i = unique(gemma_df$ee.ShortName), .options.snow = opts)
gemma_annotation_df = foreach::`%dopar%`(boot, get_gemma_annotation(i)) %>% 
  data.table::rbindlist()
parallel::stopCluster(cl)

gemma_annotation_df = gemma_annotation_df %>%
  as_tibble() %>%
  dplyr::select(ee.shortName, class.Type, class.Name, term.Name)
stopCluster(cl)
write_tsv(gemma_annotation_df, paste0('data/gemma/kegg_disease/',organism,'/annotation_summary.tsv'))
