rm(list = ls())
setwd('~/eabenchmark/') 

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(igraph))
suppressMessages(library(parallel))

option_list = list(
  make_option(c("-c", "--ncores"), type="integer", default='1',
              help=" (by default: %default)", metavar="character"),
  make_option(c("-m", "--method"), type="character", default='anubix',
              help=" (by default: %default)", metavar="character"),
  make_option(c("-n", "--network"), type="character", default='funcoup',
              help=" (by default: %default)", metavar="character"),
  make_option(c("-x", "--max_n"), type="integer", default='500',
              help=" (by default: %default)", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
method_name = opt$method
ncores = opt$ncores
network_name = opt$network
max_n = opt$max_n

source('src/others/loadNetwork.R')
source('src/others/mapGenes.R')

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

if(network_name %in% c('funcoup','string')){
  mapping_df = load_mappingFile(network_name)
  kegg_df = kegg_df %>%
    left_join(mapping_df, by = 'GeneID') %>%
    na.omit() %>%
    arrange(pathway_id) %>%
    select(Ensembl,pathway_id) %>% 
    unique()
}

################################################
### NETWORK 
graph = load_network(network_name, cutoff = 0.8)
network = igraph::as_data_frame(graph)

#################################################
## EXPRESSION DATA
dataset_df = read_tsv(file=paste0('data/TPbenchmark/cc_deg-alpha10OR20-Min15Max',max_n,'-beta0-DEG.tsv'),
                      col_types = 'fcfffiddiciccccc')

outdir = paste0('results/TPbenchmark/',max_n,'/')

### Pre-computed crosstalks
if (grepl('anubix',method_name)){
  network = network[,1:2]
  network$score = 1
  
  network_matrix = 
    ANUBIX::anubix_links(
      network,
      pathways = as.data.frame(kegg_df),
      cutoff = 0.8,
      network_type = "weighted")
  
} else if (method_name == 'binox'){
  network = network[,1:2]
  network$score = 1
  write_tsv(network, file = paste0('data/methodSpecific/binox/',network_name,'.tsv'))
  
  network_file = paste0('data/methodSpecific/binox/BinoX-',network_name,'.tsv.randNet')
  if (!file.exists(network_file)){
    ## Update path to binox
    system(paste0("/scratch/dbuzzao/PhD/binox/BinoX -n data/binox/",network_name,".tsv -s 12345 -c 0.8"))
    system(paste0("mv BinoX-",network_name,".tsv.randNet data/methodSpecific/binox/BinoX-",network_name,".tsv.randNet"))
  }
  
  out = apply(dataset_df,1,function(x){
    df = 
      data.frame(
        Ensembl = unlist(stringr::str_split(
          string = ifelse(network_name=='funcoup',x[['Ensembl']],x[['Ensembl_PRO']]),
          pattern = ',')),
        GEO = x[['GEO']]
    )
    write_tsv(df,file = paste0('data/methodSpecific/binox/',network_name,'/',x[['GEO']],'.tsv'))
    
    }) 
  rm(out)
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
## EXECUTE METHOD
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
    runtime = proc.time() - ptm
    
  } else if (method_name == 'binox'){
  ##########################################
  #### BINOX
    source('src/EAmethods/binox.R')
    # ptm = proc.time()
    tmp_enrichment_df = 
      tryCatch(
        run_binox(x[['GEO']],network_name,benchmark='positive'),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      tmp_enrichment_df = tmp_enrichment_df %>% 
          # select(`#2:NameGroupA`, `#3:NameGroupB`,`#4:p.value`,`#5:FDR`,`#14:expectedLinks`,`#15:sharedOriglinks`) %>%
          dplyr::rename(GEO = `#3:NameGroupB`,
                        pathway_id = `#2:NameGroupA`,
                        Pvalue = `#4:p.value`,
                        expected_nab = `#14:expectedLinks`,
                        nab = `#15:sharedOriglinks`) %>%
          mutate(pathway_id = tolower(pathway_id)) %>%
          select(GEO,pathway_id,nab,expected_nab,Pvalue)
  } else {return()}
    runtime = read_table('results/binox_runtime.txt', col_names=c('elapsed','user','system'), col_types='tdd') %>% 
      select(user,system,elapsed)
    
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
    
    ptm = proc.time()
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
    runtime = proc.time() - ptm
  } else if (method_name == 'netpea'){
    ##########################################
    #### netPEA
    source('src/EAmethods/netpea.R')
    Ensembl = as.character(
      unlist(stringr::str_split(
        ifelse(network_name=='funcoup',x[['Ensembl']],x[['Ensembl_PRO']]),
        ',')))

    ptm = proc.time()
    tmp_enrichment_df = 
      tryCatch(
        run_netpea(
          graph,
          Ensembl,
          kegg_df,
          1,
          x[['seed_n']]
        ),
        error=function(cond) {
          return(NULL)})
    
    if (!is.null(tmp_enrichment_df)){
      tmp_enrichment_df = tmp_enrichment_df %>% 
        add_column(GEO = x[['GEO']]) %>%
        select(GEO,pathway_id,Score,avgScore,sdScore,Z,Pvalue)
    } else {return()}
    runtime = proc.time() - ptm
    }

  tmp_enrichment_df$user=round(as.double(runtime[[1]]),2)
  tmp_enrichment_df$system=round(as.double(runtime[[2]]),2)
  tmp_enrichment_df$elapsed=round(as.double(runtime[[3]]),2) 
  ### FOR ALL METHODS
  return(tmp_enrichment_df)
}

# if(file.exists(paste0(outdir,method_name,'-',network_name,'.tsv.gz'))){
#   enrichment_df = read_tsv(paste0(outdir,method_name,'-',network_name,'.tsv.gz'), col_types = cols()) %>%
#   dataset_df = dataset_df %>%
#     filter(GEO %in% unique(enrichment_df$GEO)) %>% droplevels()
# } else {
#   enrichment_df = data.frame()
# }
# 
# print(paste0('Datasets #',nrow(dataset_df)))
# count = 1
# tot = nrow(dataset_df)

pb = utils::txtProgressBar(min=0, max=nrow(dataset_df), style = 3)
progress = function(n) utils::setTxtProgressBar(pb, n)
opts = list(progress = progress)
cl = parallel::makeCluster(ncores,outfile = paste0(method_name,'-',network_name,'.txt'))
doSNOW::registerDoSNOW(cl)
boot = foreach::foreach(i = dataset_df$GEO, .options.snow = opts)
output_df = foreach::`%dopar%`(boot, run_enrichment(i)) %>% 
  data.table::rbindlist()
parallel::stopCluster(cl)

output_df %>%
  select(-c(user,system,elapsed)) %>%
  write_tsv(file = paste0(outdir,method_name,'-',network_name,'.tsv.gz'))
# 
# output_df %>%
#   select(GEO,user,system,elapsed) %>%
#   unique() %>%
#   write_tsv(file = paste0('results/runtime/',method_name,'-',network_name,'.tsv.gz'))
