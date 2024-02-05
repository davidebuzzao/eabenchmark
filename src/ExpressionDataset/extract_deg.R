rm(list = ls())
setwd('~/eabenchmark/') 

library(tidyverse)
library(parallel)
source('src/others/loadNetwork.R')

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

funcoup_net =  load_network(network.name = 'funcoup',cutoff = 0.8)
string_net = load_network(network.name = 'string',cutoff = 0.8)

## Select TP or FP benchmark here, and run the code
benchmark = 'TPbenchmark'
# benchmark = 'FPbenchmark'

min_n = 15
max_n = 500
dataset_df = read_tsv(file=paste0('data/TPbenchmark/cc_deg-alpha10OR20-Min',min_n,'Max',max_n,'-beta0.tsv'),
                      col_types = 'fcfffiifddc')

if(benchmark=='FPbenchmark'){
  n_sampling = 30
  ## REMOVE ENSEMBL --> GENEID for >1 times
  uniprot_mapping_df = 
    uniprot_mapping_df %>% 
    filter(Ensembl %in% 
             (uniprot_mapping_df %>% 
                group_by(Ensembl) %>%
                summarise(TOT=n()) %>% 
                filter(TOT==1) %>%
                pull(Ensembl)))
  
  ## REMOVE ENSEMBL_PRO --> GENEID for >1 times
  string_mapping_df = 
    string_mapping_df %>% 
    filter(Ensembl_PRO %in% 
             (string_mapping_df %>% 
                group_by(Ensembl_PRO) %>%
                summarise(TOT=n()) %>% 
                filter(TOT==1) %>%
                pull(Ensembl_PRO)))
  
  dataset_df = dataset_df %>%
    group_by(GEO,Title,TargetPathway,DataType,BatchEffect,
             GenomeCoverage,NrOfDEG,alpha,beta,rankBy,MinMax)%>%
    summarise(
      GEO_x = list(paste(GEO,seq(n_sampling),sep = '.'))
    ) %>%
    unnest(cols = GEO_x) %>%
    ungroup() %>%
    select(-GEO) %>%
    dplyr::rename(GEO=GEO_x)
}

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

count = 1; tot = nrow(dataset_df)

extract_deg =  function(row_n){
  
  library(tidyverse)
  library(SummarizedExperiment)
  x = dataset_df[row_n,]
  se = readRDS(paste0('data/TPbenchmark/',x[['DataType']],'/',gsub('\\.[0-9]*','',x[['GEO']]),'.Rds'))
  dea_df = as.data.frame(rowData(se))
  alpha = as.double(x[['alpha']])
  beta = as.double(x[['beta']])
  num_deg = x[['NrOfDEG']]
  
  if(benchmark=='TPbenchmark'){
    # print(beta)
    if (beta>0){
      tmp_dea_df = dea_df %>% 
        filter(adj.P.Val < alpha &
                 abs(logFC) > beta) %>% 
        droplevels()
    } else {
      tmp_dea_df = dea_df %>% 
        filter(adj.P.Val < alpha) %>% 
        droplevels() %>%
        arrange(-abs(t)) %>%
        head(as.numeric(x[['NrOfDEG']]))
    }
    out_df = 
      data.frame(
        GEO = x[['GEO']],
        Title = x[['Title']], 
        TargetPathway = x[['TargetPathway']],
        DataType = x[['DataType']], 
        BatchEffect = x[['BatchEffect']],
        GenomeCoverage = x[['GenomeCoverage']],
        alpha = alpha,
        beta = beta)
    
    samples_n = table(se$GROUP)
    deg_n = table(sign(tmp_dea_df$logFC))
    out_df$NrOfSamples = as.character(sum(samples_n))
    out_df$CC = paste0(as.character(samples_n[2]),'/',as.character(samples_n[1]))
    out_df$NrOfDEG = as.character(sum(deg_n))
    out_df$UpDown = paste0(as.character(deg_n[2]),'/',as.character(deg_n[1]))
    
    out_df$GeneID = NA
    out_df$Ensembl = NA
    out_df$Ensembl_PRO = NA
    
    GeneID = unique(rownames(tmp_dea_df))

    tmp_mapping_df = mapping_df %>% 
      filter(GeneID %in% rownames(tmp_dea_df))
    
    Ensembl = unique(tmp_mapping_df$Ensembl[!is.na(tmp_mapping_df$Ensembl)])
    Ensembl_PRO = unique(tmp_mapping_df$Ensembl_PRO[!is.na(tmp_mapping_df$Ensembl_PRO)])
    
    out_df$mappedDEG = paste0(length(GeneID),'/',sum(Ensembl %in% igraph::V(funcoup_net)$name),'/',sum(Ensembl_PRO%in% igraph::V(string_net)$name))
    out_df$GeneID = paste(GeneID,collapse=",")
    out_df$Ensembl = paste(Ensembl,collapse=",")
    out_df$Ensembl_PRO = paste(Ensembl_PRO,collapse=",")
    
    return(out_df)
    
  } else if(benchmark=='FPbenchmark'){
    geo = unlist(str_split(x[['GEO']],'\\.')) 
    permuted_se = read_tsv(file = paste0('data/FPbenchmark/',geo[1],'.tsv.gz'), 
                           col_types = paste0(rep('c',30),collapse = ''))
    rownames(dea_df) = permuted_se %>% pull(paste0('P',geo[2]))
    
    if (beta>0){
      tmp_dea_df = dea_df %>% 
        filter(adj.P.Val < alpha &
                 abs(logFC) > beta) %>% 
        droplevels()
    } else {
      tmp_dea_df = dea_df %>% 
        filter(adj.P.Val < alpha) %>% 
        droplevels() %>%
        arrange(-abs(t)) %>%
        head(as.numeric(x[['NrOfDEG']]))
    }
    out_df = 
      data.frame(
        GEO = x[['GEO']],          
        Title = x[['Title']], 
        TargetPathway = x[['TargetPathway']],
        DataType = x[['DataType']], 
        BatchEffect = x[['BatchEffect']],
        GenomeCoverage = x[['GenomeCoverage']],
        alpha = alpha,
        beta = beta)
    
    samples_n = table(se$GROUP)
    deg_n = table(sign(tmp_dea_df$logFC))
    out_df$NrOfSamples = as.character(sum(samples_n))
    out_df$CC = paste0(as.character(samples_n[2]),'/',as.character(samples_n[1]))
    out_df$NrOfDEG = as.character(sum(deg_n))
    out_df$UpDown = paste0(as.character(deg_n[2]),'/',as.character(deg_n[1]))
    
    out_df$GeneID = NA
    out_df$Ensembl = NA
    out_df$Ensembl_PRO = NA
    
    GeneID = unique(rownames(tmp_dea_df))
    
    tmp_mapping_df = mapping_df %>% 
      filter(GeneID %in% rownames(tmp_dea_df))
    
    Ensembl = unique(tmp_mapping_df$Ensembl[!is.na(tmp_mapping_df$Ensembl)])
    Ensembl_PRO = unique(tmp_mapping_df$Ensembl_PRO[!is.na(tmp_mapping_df$Ensembl_PRO)])
    
    out_df$mappedDEG = paste0(length(GeneID),'/',sum(Ensembl %in% igraph::V(funcoup_net)$name),'/',sum(Ensembl_PRO%in% igraph::V(string_net)$name))
    out_df$GeneID = paste(GeneID,collapse=",")
    out_df$Ensembl = paste(Ensembl,collapse=",")
    out_df$Ensembl_PRO = paste(Ensembl_PRO,collapse=",")
    return(out_df)
  }
  
}

#################################################
## EXECUTE METHOD
pb = utils::txtProgressBar(min=0, max=nrow(dataset_df), style = 3)
progress = function(n) utils::setTxtProgressBar(pb, n)
opts = list(progress = progress)
cl = parallel::makeCluster(10)
doSNOW::registerDoSNOW(cl)
boot = foreach::foreach(i = 1:nrow(dataset_df), .options.snow = opts)
deg_df = foreach::`%dopar%`(boot, extract_deg(i)) %>% 
  data.table::rbindlist() %>%
  select(GEO, Title, TargetPathway, DataType, BatchEffect, GenomeCoverage, alpha, beta, NrOfSamples, CC, NrOfDEG, UpDown, mappedDEG, GeneID, Ensembl, Ensembl_PRO) 
parallel::stopCluster(cl)

if (benchmark == 'TPbenchmark'){
  write_tsv(deg_df %>% arrange(GEO),file = paste0('data/TPbenchmark/cc_deg-alpha10OR20-Min',min_n,'Max',max_n,'-beta0-DEG.tsv'))
  write_tsv(deg_df %>% arrange(GEO) %>% select(-c(GeneID, Ensembl, Ensembl_PRO)),file = paste0('data/TPbenchmark/cc_deg-alpha10OR20-Min',min_n,'Max',max_n,'-beta0-noDEG.tsv'))
} else {
  write_tsv(deg_df,file = paste0('data/FPbenchmark/cc_deg-alpha10OR20-Min',min_n,'Max',max_n,'-beta0-DEG.tsv.gz'))
}
