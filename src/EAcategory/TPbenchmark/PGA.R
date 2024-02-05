rm(list = ls())
setwd('~/eabenchmark/') 

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(parallel))

option_list = list(
  make_option(c("-c", "--ncores"), type="integer", default='1',
              help=" (by default: %default)", metavar="character"),
  make_option(c("-m", "--method"), type="character", default='fgsea',
              help=" (by default: %default)", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
method_name = opt$method
ncores = opt$ncores

#################################################
## FUNTIONAL GENE SETS
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
                      col_types = 'fcfffiddiciccccc') %>% 
  dplyr::rename(seed_n = seed.number)

outdir = paste0('results/TPbenchmark/')

#################################################
## EXECUTE METHOD
run_enrichment = function(d){
  library(tidyverse)
  message(d)
  #################################################
  ## EXPRESSION DATA
  x = dataset_df[dataset_df$GEO==d,] %>%
    droplevels()
  se = readRDS(paste0('data/TPbenchmark/',x[['DataType']],'/',x[['GEO']],'.Rds'))
  se@colData$GROUP = as.numeric(as.character(se@colData$GROUP))
  se@metadata$dataType = as.character(x[['DataType']])
  print(se)
  if (method_name == 'camera'){ 
   ##########################################
   #### CAMERA
   source('src/EAmethods/camera.R')
   ptm = proc.time()
   tmp_enrichment_df = 
     tryCatch(
       run_camera(se, kegg_list, var_corr = FALSE),
       error=function(cond) {
         return(NULL)})
   if (!is.null(tmp_enrichment_df)){
     p = rownames(tmp_enrichment_df)
     tmp_enrichment_df = tmp_enrichment_df %>%
       dplyr::rename(Pvalue=PValue) %>%
       magrittr::set_rownames(NULL) %>%
       add_column(Correlation=0.01) %>%
       add_column(pathway_id = p,
                  GEO = x[['GEO']]) %>%
       select(GEO,pathway_id,NGenes,Direction,Pvalue,Correlation)
   } else {return()}
   
  } else if (method_name == 'camera*'){
    ##########################################
    #### CAMERA*
    source('src/EAmethods/camera.R')
    ptm = proc.time()
    tmp_enrichment_df = 
      tryCatch(
        run_camera(se, kegg_list, var_corr = TRUE),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      p = rownames(tmp_enrichment_df)
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(Pvalue=PValue) %>%
        magrittr::set_rownames(NULL) %>%
        add_column(pathway_id = p,
                   GEO = x[['GEO']]) %>%
        select(GEO,pathway_id,NGenes,Direction,Pvalue,Correlation)
    } else {return()}
    
  } else if (method_name == 'fgsea'){
    ##########################################
    #### fGSEA
    source('src/EAmethods/fgsea.R')
    ptm = proc.time()
    tmp_enrichment_df = 
      tryCatch(
        run_fgseaMultilevel(se, kegg_list, x[['seed_n']]),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(pathway_id = pathway,
                      Pvalue = pval) %>%
        add_column(GEO = x[['GEO']]) %>%
        select(GEO,pathway_id,ES,NES,Pvalue,size,log2err,leadingEdge)
    } else {return()}
    
  } else if (method_name == 'gsea'){
    ##########################################
    #### GSEA
    source('src/EAmethods/gsea.R')
    ptm = proc.time()
    tmp_enrichment_df = 
      tryCatch(
        run_gsea(se, kegg_df, x[['seed_n']]),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      tmp_enrichment_df = tmp_enrichment_df@result %>%
        dplyr::rename(pathway_id = ID,
                      Pvalue = pvalue) %>%
        add_column(GEO = x[['GEO']]) %>%
        select(GEO,pathway_id,enrichmentScore,NES,Pvalue,setSize,leading_edge,core_enrichment)
    } else {return()}
    
  } else if (method_name == 'gsa'){
    ##########################################
    #### GSA
    source('src/EAmethods/gsa.R')
    ptm = proc.time()
    tmp_enrichment_df = 
      tryCatch(
         run_gsa(se, kegg_list, x[['seed_n']]),
         error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      tmp_enrichment_df = tmp_enrichment_df %>%
        add_column(pathway_id = names(kegg_list),
                   GEO = x[['GEO']]) %>%
        select(GEO,pathway_id,ES,Pvalue)
    } else {return()}
    
  } else if (method_name == 'gsva'){
    ##########################################
    #### GSVA
    source('src/EAmethods/gsva.R')
    ptm = proc.time()
    tmp_enrichment_df = 
      tryCatch(
        run_gsva(se, kegg_list, x[['seed_n']]),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      p = rownames(tmp_enrichment_df)
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(Pvalue = P.Value) %>%
        magrittr::set_rownames(NULL) %>%
        add_column(pathway_id = p,
                   GEO = x[['GEO']]) %>%
        select(GEO,pathway_id,t,Pvalue)
    } else {return()}
    
  } else if (method_name == 'padog'){
    ##########################################
    #### PADOG
    source('src/EAmethods/padog.R')
    ptm = proc.time()
    tmp_enrichment_df = 
      tryCatch(
        run_padog(se, kegg_list, 8, x[['seed_n']]),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      # p = rownames(tmp_enrichment_df)
      tmp_enrichment_df = tmp_enrichment_df %>%
        magrittr::set_rownames(NULL) %>%
        dplyr::rename(Pvalue = Ppadog,
                      pathway_id = ID) %>%
        add_column(GEO = x[['GEO']]) %>%
        mutate(meanAbsT0 = round(meanAbsT0,2),
               padog0 = round(padog0,4)) %>%
        select(GEO,pathway_id,meanAbsT0,padog0,PmeanAbsT,Pvalue)
    } else {return()}
    
  } else if (method_name == 'roast'){
    ##########################################
    #### ROAST
    source('src/EAmethods/roast.R')
    ptm = proc.time()
    tmp_enrichment_df = 
      tryCatch(
        run_roast(se, kegg_list, x[['seed_n']]),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      p = rownames(tmp_enrichment_df)
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(Pvalue=PValue) %>%
        magrittr::set_rownames(NULL) %>%
        add_column(pathway_id = p,
                   GEO = x[['GEO']]) %>%
        select(GEO,pathway_id,PropDown,PropUp,Direction,Pvalue)
    } else {return()}
  }
  
  runtime = proc.time() - ptm
  tmp_enrichment_df$user=round(runtime[[1]],2)
  tmp_enrichment_df$system=round(runtime[[2]],2)
  tmp_enrichment_df$elapsed=round(runtime[[3]],2)
  ### FOR ALL METHODS
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

## output runtime
output_df %>%
  select(GEO,user,system,elapsed) %>%
  unique() %>%
  write_tsv(file = paste0('results/runtime/',method_name,'.tsv.gz'))