rm(list = ls())
setwd('~/eabenchmark/') 

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(igraph))
suppressMessages(library(parallel))

source('src/others/loadNetwork.R')
source('src/others/mapGenes.R')

option_list = list(
  make_option(c("-c", "--ncores"), type="integer", default='4',
              help=" (by default: %default)", metavar="character"),
  make_option(c("-m", "--method"), type="character", default='anubix',
              help=" (by default: %default)", metavar="character"),
  make_option(c("-n", "--network"), type="character", default='funcoup',
              help=" (by default: %default)", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
method_name = opt$method
ncores = opt$ncores
network_name = opt$network

message(toupper(paste(method_name,network_name,sep = '-')))

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

if(network_name %in% c('funcoup','string')){
  mapping_df = load_mappingFile(network_name)
  kegg_df = kegg_df %>%
    left_join(mapping_df, by = 'GeneID',multiple = "all") %>%
    na.omit() %>%
    arrange(pathway_id) %>%
    select(Ensembl,pathway_id) %>% 
    unique()
}

################################################
### NETWORK 
message(paste0('Load ', toupper(network_name)))
graph = load_network(network_name, cutoff = 0.8)
network = igraph::as_data_frame(graph)

### Pre-computed crosstalks
message(paste0('Execute ', toupper(method_name)))
if (grepl('anubix',method_name)){
  network = network[,1:2]
  network$score = 1
  
  network_matrix = 
    ANUBIX::anubix_links(
      network,
      pathways = as.data.frame(kegg_df),
      cutoff = 0.8,
      network_type = "weighted")
  
} else if (method_name == 'neat'){
  network = as.matrix(network[,1:2])
  network_v = sort(igraph::V(graph)$name)
  kegg_df_filtered = kegg_df %>%
    dplyr::filter(Ensembl %in% network_v)
  kegg_list = tapply(kegg_df_filtered$Ensembl, 
                     kegg_df_filtered$pathway_id, 
                     function(x) as.character(x))
  
}

#################################################
## EXPRESSION DATA
set.seed(12345)
dataset_df = read_tsv(file='data/FPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv.gz',
                      col_types = 'fcfffiddiciccccc') %>% 
                      add_column(seed_n = sample(10000, nrow(.), replace = F))

outdir = paste0('results/FPbenchmark/')

#############################
run_enrichment = function(d){
  library(tidyverse)
  
  x = dataset_df %>% 
    filter(GEO==d) %>% droplevels()
  
  if (method_name=='anubix'){
    ##########################################
    #### ANUBIX
    source('src/EAmethods/anubix.R')
    
    gs_df = data.frame(
      Ensembl = as.character(
        unlist(stringr::str_split(
          ifelse(network_name=='funcoup',x[['Ensembl']],x[['Ensembl_PRO']]),
          ','))),
      GEO = x[['GEO']])
  
    ptm = proc.time()
    tmp_enrichment_df = 
      tryCatch(
        run_anubix(
          network, network_matrix, 
          gs_df, kegg_df,
          1, 0.8, 2000
        ),
        error=function(cond) {
          return(NULL)})

    if (!is.null(tmp_enrichment_df)){
      tmp_enrichment_df = tmp_enrichment_df %>% 
        dplyr::rename(GEO = geneset,
                      pathway_id = pathway,
                      Pvalue = `p-value`) %>%
        select(GEO,pathway_id,obv_links,exp_mean,overlap,Pvalue)
    } else {return()}
    
  } else if (method_name == 'binox'){
    ##########################################
    #### BINOX
    ## Run binox-standalone.R
    
  } else if (method_name == 'neat'){
    ##########################################
    #### NEAT
    source('src/EAmethods/neat.R')
    gs_df = data.frame(
      Ensembl = as.character(
        unlist(stringr::str_split(
          ifelse(network_name=='funcoup',x[['Ensembl']],x[['Ensembl_PRO']]),
          ','))),
      GEO = x[['GEO']])
    gs_list = stringr::str_split(ifelse(network_name=='funcoup',x[['Ensembl']],x[['Ensembl_PRO']]),',')
    names(gs_list) = x[['GEO']]
    tmp_enrichment_df = 
      tryCatch(
        as.data.frame(
          run_neat(network, network_v,
                   gs_list, kegg_list)
        ),
        error=function(cond) {
          return(NULL)})
    
    if (!is.null(tmp_enrichment_df)){
      tmp_enrichment_df = tmp_enrichment_df %>% 
        dplyr::rename(GEO = A,
                      pathway_id = B,
                      Pvalue = pvalue) %>%
        select(GEO,pathway_id,nab,expected_nab,Pvalue)
    } else {return()}
    
  } else if (method_name == 'netpea'){
    ##########################################
    #### netPEA
    source('src/EAmethods/netpea.R')
    GeneID = as.character(unlist(stringr::str_split(x[['GeneID']],',')))
    
    tmp_enrichment_df = 
      tryCatch(
        run_netpea(
          graph,
          GeneID,
          kegg_df,
          mapping_df,
          x[['seed_n']]
        ),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      tmp_enrichment_df = tmp_enrichment_df %>% 
        add_column(GEO = x[['GEO']]) %>%
        select(GEO,pathway_id,Score,avgScore,sdScore,Z,Pvalue)
    } else {return()}
  }
  
  ### FOR ALL METHODS
  return(tmp_enrichment_df)
}

print(paste0('Datasets #',nrow(dataset_df)))

pb = utils::txtProgressBar(min=0, max=nrow(dataset_df), style = 3)
progress = function(n) utils::setTxtProgressBar(pb, n)
opts = list(progress = progress)
cl = parallel::makeCluster(ncores,outfile = paste0(method_name,'-',network_name,'.txt'))
doSNOW::registerDoSNOW(cl)
boot = foreach::foreach(i = dataset_df$GEO, .options.snow = opts)
enrichment_df = foreach::`%dopar%`(boot, run_enrichment(i)) %>% 
  data.table::rbindlist()
parallel::stopCluster(cl)

enrichment_df %>% 
  arrange(GEO) %>%
  write_tsv(file = paste0(outdir,method_name,'-',network_name,'.tsv.gz'))