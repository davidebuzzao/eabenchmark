rm(list = ls())
setwd('~/eabenchmark/') 

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(parallel))

option_list = list(
  make_option(c("-c", "--ncores"), type="integer", default='1',
              help=" (by default: %default)", metavar="character"),
  make_option(c("-m", "--method"), type="character", default='spia',
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

### Pre-computed KEGG networks
if (method_name == 'spia'){
  library(SPIA)
  ## parsing up-to-date KEGG xml files for use with SPIA
  ## Download xml KEGG file 
  for (p in names(kegg_list)){
    if(file.exists(paste0('data/KEGG/xml/',p,'.xml'))) next
    system(paste0('bash -c \'',
                  'curl -s http://rest.kegg.jp/get/',p,'/kgml -o data/KEGG/xml/',p,'.xml\''))
  }
  ## or if you don't know the pathway names
  # curl -s http://rest.kegg.jp/list/pathway/hsa | awk '{split($1,a,':'); print 'curl -s http://rest.kegg.jp/get/'a[2]'/kgml -o data/kegg/xml/'a[2]'.xml'}' | bash
  if(!file.exists('data/methodSpecific/spia/hsaSPIA.RData')){
    ## Build SPIA data takes 56'
    makeSPIAdata(kgml.path='data/KEGG/xml',organism='hsa',out.path='data/methodSpecific/spia')
  }

  relation_types = 
    c('activation', 'compound', 'binding/association', 
      'expression', 'inhibition', 'activation_phosphorylation', 
      'phosphorylation', 'inhibition_phosphorylation', 
      'inhibition_dephosphorylation', 'dissociation', 'dephosphorylation', 
      'activation_dephosphorylation', 'state change', 'activation_indirect effect', 
      'inhibition_ubiquination', 'ubiquination', 'expression_indirect effect', 
      'inhibition_indirect effect', 'repression', 'dissociation_phosphorylation', 
      'indirect effect_phosphorylation', 'activation_binding/association', 
      'indirect effect', 'activation_compound', 'activation_ubiquination')
  beta_v = c(1, 0, 0, 1, -1, 1, 0, -1, -1, 0, 0, 1, 0, 1, -1, 0, 1, -1, -1, 0, 0, 1, 0, 1, 1)
  names(beta_v) = relation_types
  print(cbind(beta_v))
  
} else if (method_name %in% c('cepaORA')){
  library(CePa)
  kegg_net = read_tsv('data/KEGG/hsaKEGG_network.tsv', col_types = 'ccffffc') %>%
    dplyr::rename(input=from,
                  output=to) %>%
    group_by(input,output) %>%
    mutate(interaction.id = as.character(cur_group_id())) 
  
  kegg_list = tapply(kegg_net$interaction.id, kegg_net$pathway_id, function(x) as.character(x))
  
  kegg_net = kegg_net %>%
    select(interaction.id,input,output) %>%
    unique()
  
  mapping_df = data.frame(node.id = unique(c(kegg_net$input,kegg_net$output)),
                          GeneID = unique(c(kegg_net$input,kegg_net$output)))
  
  kegg_catalogue = 
    set.pathway.catalogue(
      pathList = kegg_list,
      interactionList = as.data.frame(kegg_net),
      mapping = mapping_df
    )
}


#################################################
## EXPRESSION DATA
dataset_df = read_tsv(file='data/TPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv', 
                      col_types = 'fcfffiddiciccccc')

outdir = paste0('results/TPbenchmark/')

#################################################
## EXECUTE METHOD
run_enrichment = function(d){
    library(tidyverse)
    x = dataset_df[dataset_df$GEO==d,] %>%
      droplevels()
    se = readRDS(paste0('data/TPbenchmark/',x[['DataType']],'/',x[['GEO']],'.Rds'))
    se@colData$GROUP = as.numeric(as.character(se@colData$GROUP))
    se@metadata$dataType = as.character(x[['DataType']])
    GeneID = as.character(unlist(stringr::str_split(x[['GeneID']],',')))
    
    if (method_name == 'spia'){
      ##########################################
      #### SPIA
      source('src/EAmethods/spia.R')
      ptm = proc.time()
      tmp_enrichment_df = 
        tryCatch(
          run_spia(se,GeneID,beta_v,x[['seed_n']]),
          error=function(cond) {
            return(NULL)})
      if (!is.null(tmp_enrichment_df)){
        tmp_enrichment_df = tmp_enrichment_df %>%
          mutate(ID = paste0('hsa',ID)) %>%
          dplyr::rename(Pvalue=pG,
                        pathway_id=ID) %>%
          add_column(GEO = x[['GEO']]) %>%
          select(GEO,pathway_id,pSize,pNDE,tA,pPERT,Pvalue,Status,KEGGLINK)
        
      } else {return()}
    } else if (method_name == 'cepaORA'){
      ##########################################
      #### CePa ORA
      source('src/EAmethods/cepa.R')
      ptm = proc.time()
      tmp_enrichment_df = 
        tryCatch(
          run_cepa(se, GeneID, kegg_catalogue, method = 'ora', x[['seed_n']]),
          error=function(cond) {
            return(NULL)})
      if (!is.null(tmp_enrichment_df)){
        kegg_names = rownames(p.table(tmp_enrichment_df))
        tmp_enrichment_df = as.data.frame(p.table(tmp_enrichment_df)) %>%
          add_column(pathway_id = kegg_names,
                     GEO = x[['GEO']]) %>%
          magrittr::set_rownames(NULL) %>%
          select(GEO,pathway_id,colnames(.))
      } else {return()}
      
    } else if (method_name == 'cepaGSA'){
      ##########################################
      #### CePa GSA
      source('src/EAmethods/cepa.R')
      ptm = proc.time()
      tmp_enrichment_df = 
        tryCatch(
          run_cepa(se, NULL, kegg_catalogue, method = 'gsa', x[['seed_n']]),
          error=function(cond) {
            return(NULL)})
      if (!is.null(tmp_enrichment_df)){
        kegg_names = rownames(p.table(tmp_enrichment_df))
        tmp_enrichment_df = as.data.frame(p.table(tmp_enrichment_df)) %>%
          add_column(pathway_id = kegg_names,
                     GEO = x[['GEO']]) %>%
          magrittr::set_rownames(NULL) %>%
          select(GEO,pathway_id,colnames(.))
      } else {return()}
    }
    
    runtime = proc.time() - ptm
    tmp_enrichment_df$user=round(runtime[[1]],2)
    tmp_enrichment_df$system=round(runtime[[2]],2)
    tmp_enrichment_df$elapsed=round(runtime[[3]],2)
    return(tmp_enrichment_df)
}

## Proceed with computation of shortest path
## Set up parallelization, either on 1 or more cores
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
