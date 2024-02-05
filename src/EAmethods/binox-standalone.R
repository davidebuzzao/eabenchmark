rm(list = ls())
indir = '/scratch2/dbuzzao/eabenchmark/'
setwd(indir) 

############################################
## 183.5' to randomize FunCoup 2000times
## 73.6' to randomize STRING 2000times
## Rscript local/src/methods/binox-standalone.R -n funcoup  -b TPbenchmark
## Rscript local/src/methods/binox-standalone.R -n funcoup -b FPbenchmark
## Rscript local/src/methods/binox-standalone.R -n string -b TPbenchmark
## Rscript local/src/methods/binox-standalone.R -n string -b FPbenchmark

suppressMessages(library(optparse))
option_list = list(
  make_option(c("-n", "--network"), type="character", default='funcoup',
              help=" (by default: %default)", metavar="character"),
  make_option(c("-b", "--benchmark"), type="character", default='TPbenchmark',
              help=" (by default: %default)", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
network_name = opt$network
benchmark = opt$benchmark

suppressMessages(library(tidyverse))
suppressMessages(library(igraph))
source('src/others/loadNetwork.R')
source('src/others/mapGenes.R')

outdir = paste0('results/',benchmark,'/')
binox_dir = '/scratch2/dbuzzao/binox/'

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

if(network_name %in% c('funcoup','string')){
  mapping_df = load_mappingFile(network_name)
  kegg_df %>%
    left_join(mapping_df, by = 'GeneID') %>%
    na.omit() %>%
    arrange(pathway_id) %>%
    select(Ensembl,pathway_id) %>%
    unique() %>% 
    write_tsv(paste0('data/KEGG/kegg-',network_name,'.tsv'))
}

################################################
### NETWORK 
network_file = paste0('data/methodSpecific/binox/BinoX-',network_name,'.tsv.randNet')
if (!file.exists(network_file)){
  graph = load_network(network_name, cutoff = 0.8)
  network = igraph::as_data_frame(graph)
  network = network[,1:2]
  network$score = 1
  write_tsv(network, file = paste0('data/methodSpecific/binox/',network_name,'.tsv'))
  rm(network)
  
  system(paste0(binox_dir,"BinoX -n ", indir, "data/methodSpecific/binox/",network_name,".tsv -s 12345 -c 0.8 -i 2000"))
  system(paste0("mv BinoX-",network_name,".tsv.randNet ", indir, "data/methodSpecific/binox/BinoX-",network_name,".tsv.randNet"))
}

if (benchmark=='TPbenchmark'){
  dataset_df = read_tsv(file='data/TPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv',
                        col_types = 'fcfffiddiciccccc') 
} else if (benchmark=='FPbenchmark') {
  # benchmark = 'negative-newAssignment'
  benchmark = 'negative-resampling'
  dataset_df = read_tsv(file='data/FPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv.gz',
                        col_types = 'fcfffiddiciccccc')
}

dataset_df = apply(dataset_df,1,function(x){
  df =
    data.frame(
      Ensembl = unlist(stringr::str_split(
        string = ifelse(network_name=='funcoup',x[['Ensembl']],x[['Ensembl_PRO']]),
        pattern = ',')),
      GEO = x[['GEO']]
    )
}) %>% data.table::rbindlist()
write_tsv(dataset_df,file = paste0('data/methodSpecific/binox/',benchmark,'_gs-',network_name,'.tsv'))

#################################################
## EXECUTE METHOD
system(paste0('bash -c \'',
              binox_dir, 'BinoX -r ',network_file,' ',
              '-a <(tail -n +2 ', indir, 'data/methodSpecific/binox/kegg-',network_name,'.tsv | cut -f 1,2) ',
              '-b <(tail -n +2 ', indir, 'data/methodSpecific/binox/',benchmark,'_gs-',network_name,'.tsv | cut -f 1,2) ',
              '-o ', indir, 'results/',benchmark,'_binox-',network_name,'.tsv ',
              '-p large ',
              '-g 1\''))

tmp_enrichment_df = read_tsv(file = paste0('results/',benchmark,'_binox-',network_name,'.tsv'))
enrichment_df = tmp_enrichment_df %>%
  dplyr::rename(GEO = `#3:NameGroupB`,
                pathway_id = `#2:NameGroupA`,
                Pvalue = `#4:p.value`,
                expected_nab = `#14:expectedLinks`,
                nab = `#15:sharedOriglinks`) %>%
  mutate(pathway_id = tolower(pathway_id)) %>%
  select(GEO,pathway_id,nab,expected_nab,Pvalue)

write_tsv(enrichment_df, paste0(outdir,'binox-',network_name,'.tsv.gz'))
