rm(list = ls())
setwd('~/eabenchmark/') 

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(parallel))

option_list = list(
  make_option(c("-c", "--ncores"), type="integer", default='8',
              help=" (by default: %default)", metavar="character"),
  make_option(c("-m", "--method"), type="character", default='ease',
              help=" (by default: %default)", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
method_name = opt$method
ncores = opt$ncores

message(paste0('Execute ', toupper(method_name)))

#################################################
## FUNTIONAL GENE SETS
message('Load KEGG')
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
  select(pathway_id,entrezgene_id) %>%
  dplyr::rename(GeneID=entrezgene_id)

kegg_list = tapply(kegg_df$GeneID, kegg_df$pathway_id, function(x) as.character(x))
kegg_list = lapply(kegg_list, function(x) unique(x))

#################################################
## EXPRESSION DATA
set.seed(12345)
dataset_df = read_tsv(file='data/FPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv.gz',
                      col_types = 'fcfffiddiciccccc') %>%
  add_column(seed_n = sample(10000, nrow(.), replace = F))

outdir = paste0('results/FPbenchmark/')

#################################################
## EXECUTE METHOD
run_enrichment = function(d){
  library(tidyverse)
  x = dataset_df[dataset_df$GEO==d,] %>%
    droplevels()
  
  if (method_name == 'ease'){
    source('src/EAmethods/ease.R')
    tmp_enrichment_df = run_ease(as.character(unlist(stringr::str_split(x[['GeneID']],','))),
                                 PT=as.numeric(x[['GenomeCoverage']]),
                                 kegg_list)
    tmp_enrichment_df$GEO = x[['GEO']]
    return(tmp_enrichment_df)
  } else if (method_name == 'fisher'){
    source('src/EAmethods/fisher.R')
    tmp_enrichment_df = run_fisher(as.character(unlist(stringr::str_split(x[['GeneID']],','))),
                                   PT=as.numeric(x[['GenomeCoverage']]),
                                   kegg_list)
    tmp_enrichment_df$GEO = x[['GEO']]
    return(tmp_enrichment_df)
  } else if (method_name == 'random'){
    source('src/EAmethods/random.R')
    tmp_enrichment_df = run_random(as.character(unlist(stringr::str_split(x[['GeneID']],','))),
                                   kegg_list)
    tmp_enrichment_df$GEO = x[['GEO']]
  }
  
  ### FOR ALL METHODS
  return(tmp_enrichment_df)
}

pb = utils::txtProgressBar(min=0, max=nrow(dataset_df), style = 3)
progress = function(n) utils::setTxtProgressBar(pb, n)
opts = list(progress = progress)
cl = parallel::makeCluster(ncores)
doSNOW::registerDoSNOW(cl)
boot = foreach::foreach(i = dataset_df$GEO, .options.snow = opts)
enrichment_df = foreach::`%dopar%`(boot, run_enrichment(i)) %>% 
  data.table::rbindlist()
parallel::stopCluster(cl)

enrichment_df %>%
  write_tsv(file = paste0(outdir,method_name,'.tsv.gz'))
