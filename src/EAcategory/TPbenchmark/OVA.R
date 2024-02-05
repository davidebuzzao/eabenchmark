rm(list = ls())
setwd('~/eabenchmark/') 

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(parallel))

option_list = list(
  make_option(c("-c", "--ncores"), type="integer", default='1',
              help=" (by default: %default)", metavar="character"),
  make_option(c("-m", "--method"), type="character", default='ease',
              help=" (by default: %default)", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
method_name = opt$method
ncores = opt$ncores

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
dataset_df = read_tsv(file='data/TPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv', 
                      col_types = 'fcfffiddiciccccc')

outdir = paste0('results/TPbenchmark/')

#################################################
## EXECUTE METHOD
run_enrichment = function(d){
  library(tidyverse)
  
  x = dataset_df[dataset_df$GEO==d,]
  
  if (method_name == 'ease'){
    source('src/EAmethods/ease.R')
    ptm = proc.time()
    tmp_enrichment_df = run_ease(as.character(unlist(stringr::str_split(x[['GeneID']],','))),
                                 PT=as.numeric(x[['GenomeCoverage']]),
                                 kegg_list)
    tmp_enrichment_df$GEO = x[['GEO']]
  } else if (method_name == 'fisher'){
    source('src/EAmethods/fisher.R')
    ptm = proc.time()
    tmp_enrichment_df = run_fisher(as.character(unlist(stringr::str_split(x[['GeneID']],','))),
                                   PT=as.numeric(x[['GenomeCoverage']]),
                                   kegg_list)
    tmp_enrichment_df$GEO = x[['GEO']]
  } else if (method_name == 'random'){
    source('src/EAmethods/random.R')
    ptm = proc.time()
    tmp_enrichment_df = run_random(as.character(unlist(stringr::str_split(x[['GeneID']],','))),
                                   kegg_list)
    tmp_enrichment_df$GEO = x[['GEO']]
  }
  
  ### FOR ALL METHODS
  runtime = proc.time() - ptm
  tmp_enrichment_df$user=round(runtime[[1]],2)
  tmp_enrichment_df$system=round(runtime[[2]],2)
  tmp_enrichment_df$elapsed=round(runtime[[3]],2)
  return(tmp_enrichment_df)
}

pb = utils::txtProgressBar(min=0, max=nrow(dataset_df), style = 3)
progress = function(n) utils::setTxtProgressBar(pb, n)
opts = list(progress = progress)
cl = parallel::makeCluster(ncores)
doSNOW::registerDoSNOW(cl)
boot = foreach::foreach(i = dataset_df$GEO, .options.snow = opts)
output_df = foreach::`%dopar%`(boot, run_enrichment(i)) %>% 
  data.table::rbindlist()
parallel::stopCluster(cl)

## output TPbenchmark
output_df %>%
  select(-c(user,system,elapsed)) %>%
  write_tsv(file = paste0(outdir,method_name,'.tsv.gz'))

# output runtime
output_df %>%
  select(GEO,user,system,elapsed) %>%
  unique() %>%
  write_tsv(file = paste0('results/runtime/',method_name,'.tsv.gz'))
