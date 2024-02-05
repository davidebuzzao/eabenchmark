#!/usr/bin/env Rscript

rm(list = ls())
setwd('~/eabenchmark/') 

suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(parallel))
source('src/others/loadNetwork.R')

option_list = list(
  make_option(c("-c", "--ncores"), type="integer", default='4',
              help=" (by default: %default)", metavar="character"),
  make_option(c("-n", "--network"), type="character", default='humanet',
              help=" (by default: %default)", metavar="character"),
  make_option(c("-k", "--chunk"), type="integer", default=1,
              help=" (by default: %default)", metavar="integer")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
ncores = opt$ncores
network_name = opt$network
chunk = opt$chunk

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
net_edges = igraph::get.edgelist(net_graph)

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

network_degree %>% 
  group_by(bin) %>%
  count() %>% 

d = network_degree %>% 
  group_by(bin) %>% 
  summarise(diff = max(degree)-min(degree))

connectivity_df = 
  function(i){
    library(tidyverse)
    extract_net = 
      function(source_v){
        subnet_edges = net_edges[(net_edges[,1] %in% source_v | net_edges[,2] %in% source_v),]
        unique(c(subnet_edges[,1],subnet_edges[,2]))
      }
    
    compute_ji = 
      function(A,B){
        return(length(intersect(A,B))/length(union(A,B)))
      }
    p = pathway_combinations[i,]
    
    # extract links between the vertex set and the rest of the network, 
    # excluding links between the vertices in the set
    
    A_name = p[['P1']]
    B_name = p[['P2']]
    
    A_genes = kegg_list[[A_name]]
    B_genes = kegg_list[[B_name]]
    A_size = length(A_genes)
    B_size = length(B_genes)
    
    A = extract_net(A_genes)
    ## Exclude links within A and to B
    A = A[!A %in% setdiff(A_genes,intersect(A_genes,B_genes))]
    
    B = extract_net(B_genes)
    ## Exclude links within B and to A
    B = B[!B %in% setdiff(B_genes,intersect(B_genes,A_genes))]
    
    ## Observed Xtalk
    Xtalk = compute_ji(A,B)
    
    if(is.nan(Xtalk)){
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
    B = extract_net(B_genes)
    sample_XtalkAB =
      sapply(1:n_samples, function(i){
        S = sample_A[sample_A$sample_num==i,]
        S_genes = S$sample_GeneID
        S = extract_net(S_genes)
        ## Exclude links within S and to B
        S = S[!S %in% setdiff(S_genes,intersect(S_genes,B_genes))]
        ## Exclude links within B and to S
        tmp_B = B[!B %in% setdiff(S_genes,intersect(S_genes,B_genes))]
        ## Observed shared Xtalk
        compute_ji(S,tmp_B)
      })
    
    # permutation test
    sample_pvalAB = (sum(sample_XtalkAB>=Xtalk)+1)/(n_samples+1)
    ## Assume normal distribution
    expected_Xtalk_AB = mean(sample_XtalkAB)
    sd_Xtalk_AB = sd(sample_XtalkAB)
    zscore_AB = (Xtalk - expected_Xtalk_AB) / sd_Xtalk_AB
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
    A = extract_net(A_genes)
    sample_XtalkBA =
      sapply(1:n_samples, function(i){
        S = sample_A[sample_A$sample_num==i,]
        S_genes = S$sample_GeneID
        S = extract_net(S_genes)
        ## Exclude links within S and to A
        S = S[!S %in% setdiff(S_genes,intersect(S_genes,A_genes))]
        ## Exclude links within A and to S
        tmp_A = A[!A %in% setdiff(S_genes,intersect(S_genes,A_genes))]
        ## Observed shared Xtalk
        compute_ji(S,tmp_A)
      })
    
    # permutation test
    sample_pvalBA = (sum(sample_XtalkBA>=Xtalk)+1)/(n_samples+1)
    ## Assume normal distribution
    expected_Xtalk_BA = mean(sample_XtalkBA)
    sd_Xtalk_BA = sd(sample_XtalkBA)
    zscore_BA = (Xtalk - expected_Xtalk_BA) / sd_Xtalk_BA
    norm_pvalBA = pnorm(zscore_BA, mean = 0, sd = 1, lower.tail = FALSE)
    
    return(
      data.frame(
        ##data
        A = A_name,
        B = B_name,
        direction = c('AB','BA'),
        ##observation
        Xtalk = Xtalk,
        sample_avg_Xtalk = c(expected_Xtalk_AB,expected_Xtalk_BA),
        sample_sd_Xtalk = c(sd_Xtalk_AB,sd_Xtalk_BA),
        ##test
        sample_pval = c(sample_pvalAB,sample_pvalBA),
        zscore = c(zscore_AB,zscore_BA),
        norm_pval = c(norm_pvalAB,norm_pvalBA))
    )
  } 

message('COMPUTING Inter Pahtway Connectivity for #',nrow(pathway_combinations),' of genesets')

cl = parallel::makeCluster(ncores)
doSNOW::registerDoSNOW(cl)
pb = utils::txtProgressBar(min=0, max=nrow(pathway_combinations), style = 3)
progress = function(n) utils::setTxtProgressBar(pb, n)
opts = list(progress = progress)
boot = foreach::foreach(i = seq(nrow(pathway_combinations)), .options.snow = opts)
connectivity_df = foreach::`%dopar%`(boot, test_interpathwayconnectivity(i)) %>% 
  data.table::rbindlist()
parallel::stopCluster(cl)

write_tsv(connectivity_df, file = 'results/KEGG/InterPahtwayConnectivity_degree-aware1kSampling_KEGGinHumanNet20p.tsv.gz')
