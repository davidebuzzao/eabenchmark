#!/usr/bin/env Rscript

rm(list = ls())
setwd('~/eabenchmark/') 

suppressMessages(library(tidyverse))
suppressMessages(library(igraph))
suppressMessages(library(parallel))
source('src/others/loadNetwork.R')

option_list = list(
  make_option(c("-c", "--ncores"), type="integer", default='4',
              help=" (by default: %default)", metavar="character"),
  make_option(c("-n", "--network"), type="character", default='humanet',
              help=" (by default: %default)", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
method_name = opt$method
ncores = opt$ncores
network_name = opt$network

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
  dplyr::rename(GeneID=entrezgene_id)

kegg_list = tapply(kegg_df$GeneID, kegg_df$pathway_id, function(x) as.character(x))
kegg_list = lapply(kegg_list, function(x) unique(x))

pathway_disease = kegg_df %>%
  select(pathway_id,class,subclass) %>%
  filter(class=='Human Diseases') %>%
  unique() %>%
  pull(pathway_id) %>% 
  as.character()

## sample from KEGG genes
pool_genes = unique(unlist(kegg_list, use.names = FALSE))

## create combinations
message('EXTRACTING pathway combinations')
tmp_pathway_combinations = 
  # 318*317/2
  lapply(
    utils::combn(unique(as.character(kegg_df$pathway_id)), 2, simplify = F),
    function(x){
      return(data.frame(P1 = x[[1]],
                        P2 = x[[2]]))
    }) %>% data.table::rbindlist()

disease_only = T
if (disease_only){
  pathway_combinations = tmp_pathway_combinations %>% 
    filter(P1 %in% pathway_disease | P2 %in% pathway_disease)
}

n_samples = 1000
set.seed(12345)
pathway_combinations = pathway_combinations %>%
  add_column(seed_n = as.numeric(sample(100000,nrow(.),replace = F)))

#################################################
## NETWORK
min_bin_size = 100 # For degree binning
net_graph = load_network(network_name,cutoff='20p')

bins = igraph::degree(net_graph) %>%
  enframe(name = "node", value = "degree") %>%
  filter(node %in% pool_genes) %>%
  arrange(degree, as.numeric(node)) %>%
  mutate(bin = ceiling(row_number()/min_bin_size))


if(network_name %in% c('funcoup','string')){
  colnames(bins) = c('Ensembl','degree','bin')
  network_degree = read_tsv(paste0('data/network/',network_name,'_Ensembl-GeneID.tsv'), col_types = 'cc')
  network_degree = network_degree %>% 
    left_join(bins, by = 'Ensembl') %>%
    arrange(bin)
} else {
  colnames(bins) = c('GeneID','degree','bin')
  network_degree = bins %>% arrange(bin)
}

test_overlap = 
  function(i){
    library(tidyverse)
    
    p = pathway_combinations[i,]
    compute_ji = 
      function(A,B){
        return(length(intersect(A,B))/length(union(A,B)))
      }
    
    A_name = p[['P1']]
    A_genes = kegg_list[[A_name]]
    A_size = length(A_genes)
    
    B_name = p[['P2']]
    B_genes = kegg_list[[B_name]]
    B_size = length(B_genes)
    
    f = paste0('data/KEGG/overlap/',A_name,'-',B_name,'.tsv')
    if(file.exists(f)){return()}
    
    ## Observed overlap
    overlap = compute_ji(A_genes,B_genes)
    
    if(is.nan(overlap)){
      return()
    }
    
    ########################
    ## A --> B
    set.seed(p[['seed_n']])
    print('A-->B')
    seeds_numbers = sample(100000,A_size)
    names(seeds_numbers) = A_genes
    
    ## prepare random samples
    sample_A =
      lapply(A_genes,function(g){
        set.seed(seeds_numbers[g])
        
        if(!g %in% network_degree$GeneID){g=network_degree$GeneID[[1]]}
        bin_num = network_degree %>% 
          filter(GeneID==g) %>%
          select(bin)
        
        network_degree %>%
          filter(bin == bin_num[[1]]) %>%
          filter(GeneID != g) %>%
          sample_n(size = n_samples, replace = TRUE) %>%
          select(GeneID) %>%
          rename(sample_GeneID = GeneID) %>%
          add_column(sample_num = 1:n_samples) %>%
          select(sample_num,sample_GeneID)
        
      }) %>% data.table::rbindlist()
    
    ## random
    sample_overlapAB =
      sapply(1:n_samples, function(i){
        S = sample_A[sample_A$sample_num==i,]
        compute_ji(S$sample_GeneID,B_genes)
      })

    # permutation test
    sample_pvalAB = (sum(sample_overlapAB>=overlap)+1)/(n_samples+1)
    ## Assume normal distribution
    expected_overlap_AB = mean(sample_overlapAB)
    sd_overlap_AB = sd(sample_overlapAB)
    zscore_AB = (overlap - expected_overlap_AB) / sd_overlap_AB
    norm_pvalAB = pnorm(zscore_AB, mean = 0, sd = 1, lower.tail = FALSE)
    
    ########################
    ## B --> A
    set.seed(p[['seed_n']])
    print('B-->A')
    
    seeds_numbers = sample(100000,B_size)
    names(seeds_numbers) = B_genes
    
    ## prepare random samples
    sample_B =
      lapply(B_genes,function(g){
        set.seed(seeds_numbers[g])
        
        if(!g %in% network_degree$GeneID){g=network_degree$GeneID[[1]]}
        bin_num = network_degree %>% 
          filter(GeneID==g) %>%
          select(bin)
        
        network_degree %>%
          filter(bin == bin_num[[1]]) %>%
          filter(GeneID != g) %>%
          sample_n(size = n_samples, replace = TRUE) %>%
          select(GeneID) %>%
          rename(sample_GeneID = GeneID) %>%
          add_column(sample_num = 1:n_samples) %>%
          select(sample_num,sample_GeneID)
      }) %>% data.table::rbindlist()
    
    ## random
    sample_overlapBA =
      sapply(1:n_samples, function(i){
        S = sample_B[sample_B$sample_num==i,]
        compute_ji(S$sample_GeneID,A_genes)
      })
    
    # permutation test
    sample_pvalBA = (sum(sample_overlapBA>=overlap)+1)/(n_samples+1)
    ## Assume normal distribution
    expected_overlap_BA = mean(sample_overlapBA)
    sd_overlap_BA = sd(sample_overlapBA)
    zscore_BA = (overlap - expected_overlap_BA) / sd_overlap_BA
    norm_pvalBA = pnorm(zscore_BA, mean = 0, sd = 1, lower.tail = FALSE)
    
    return(
      data.frame(
      ##data
      A = A_name,
      B = B_name,
      direction = c('AB','BA'),
      ##observation
      overlap = overlap,
      sample_avg_overlap = c(expected_overlap_AB,expected_overlap_BA),
      sample_sd_overlap = c(sd_overlap_AB,sd_overlap_BA),
      ##test
      sample_pval = c(sample_pvalAB,sample_pvalBA),
      zscore = c(zscore_AB,zscore_BA),
      norm_pval = c(norm_pvalAB,norm_pvalBA))
    )
  } 

message('COMPUTING Inter Pathway Overlap for #',nrow(pathway_combinations),' of genesets')

pb = utils::txtProgressBar(min=0, max=nrow(pathway_combinations), style = 3)
progress = function(n) utils::setTxtProgressBar(pb, n)
opts = list(progress = progress)
cl = parallel::makeCluster(10,outfile = 'overlap_kegg.txt')
doSNOW::registerDoSNOW(cl)
boot = foreach::foreach(i = seq(nrow(pathway_combinations)), .options.snow = opts)
overlap_df = foreach::`%dopar%`(boot, test_overlap(i)) %>%
  data.table::rbindlist()
parallel::stopCluster(cl)

write_tsv(overlap_df, file = paste0('results/KEGG/InterPahtwayOverlap_degree-aware1kSampling_KEGGinHumanNet20p.tsv.gz'))

