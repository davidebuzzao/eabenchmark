## Follow instructions as in https://yulab-smu.top/biomedical-knowledge-mining-book/semantic-similarity-overview.html
## GO semantic similarity is an expensive operation. To be run on computer cluster. 
# Rscript local/src/diseasenetwork/semantic-genesets_analysis.R --ontology BP --ncores 120
rm(list = ls())
setwd('~/eabenchmark/') 

library(tidyverse)
library(parallel)
library(optparse)

## To use from RStudio
option_list = list(
  make_option(c("-c", "--ncores"), type="integer", default='120',
              help=" (by default: %default)", metavar="character"),
  make_option(c("-o", "--ontology"), type="character", default='BP',
              help=" (by default: %default)", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
ncores = opt$ncores
ontology = opt$ontology

fp = paste0('data/KEGG/semantic/',ontology)
if (!dir.exists(fp)) {
  dir.create(fp,recursive = T)
  print("Folder created.")
} else {
  print("Folder already exists.")
}

#########################################
########## SECTION 1 - Semantic analysis
# Functional similarity of gene products can be estimated by controlled biological 
# vocabularies, such as Gene Ontology (GO). Here we compute GO semantic similarity 
# only for significant pathway-pathway relationships

message(paste0('## COMPUTE GO-', ontology,' SEMANTIC SIMILARITY'))

## Inter-pathway overlap
overlap_df = 
  read_tsv('results/KEGG/InterPahtwayOverlap_degree-aware1kSampling_KEGGinHumanNet20p.tsv.gz', 
           col_select = c('A','B','direction','overlap','sample_pval'), 
           col_types = 'ccfdd') %>%
  dplyr::rename(overlap_pval = sample_pval) %>%
  # Swap and sort columns
  rowwise() %>%
  mutate(A = sort(c(A, B))[1],
         B = sort(c(A, B))[2]) %>%
  ungroup()

## Inter-pathway connectivity
xtalk_df = read_tsv('results/KEGG/InterPahtwayConnectivity_degree-aware1kSampling_KEGGinHumanNet20p.tsv.gz', 
                    col_select = c('A','B','direction','Xtalk','sample_pval'), 
                    col_types = 'ccfdd') %>%
  dplyr::rename(xtalk_pval = sample_pval) %>%
  # Swap and sort columns
  rowwise() %>%
  mutate(A = sort(c(A, B))[1],
         B = sort(c(A, B))[2]) %>%
  ungroup()

### TO compare 
overlap_xtalk_df = overlap_df %>% 
  inner_join(xtalk_df, by = c('A','B','direction'), multiple = 'all') %>% 
  group_by(A,B) %>% 
  summarise(
    overlap_sumlog = metap::sumlog(overlap_pval)[['p']],
    xtalk_sumlog = metap::sumlog(xtalk_pval)[['p']]) %>%
  ungroup() %>%
  mutate(
    overlap_FDR = p.adjust(overlap_sumlog, method = 'BH'),
    xtalk_FDR = p.adjust(xtalk_sumlog, method = 'BH')) %>%
  filter(overlap_FDR<0.05 & xtalk_FDR<0.05) %>% 
  select(A,B)


compute_GOsemantic = 
  function(overlap_xtalk_df_chunk){
    library(tidyverse)
    library(GOSemSim)

    if (ontology=='BP'){
      hsGO = godata('org.Hs.eg.db', ont="BP", computeIC=F)
    } else if (ontology=='MF'){
      hsGO = godata('org.Hs.eg.db', ont="MF", computeIC=F)
    } else if (ontology=='CC'){
      hsGO = godata('org.Hs.eg.db', ont="CC", computeIC=F)
    }
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
      dplyr::rename(GeneID=entrezgene_id)
    
    kegg_list = tapply(kegg_df$GeneID, kegg_df$pathway_id, function(x) as.character(x))
    kegg_list = lapply(kegg_list, function(x) unique(x))
    
    for (i in 1:nrow(overlap_xtalk_df_chunk)){
      x = overlap_xtalk_df_chunk[i,]
      f = paste0('data/KEGG/semantic/',ontology,'/',x[['A']],'-',x[['B']],'.tsv')
      if(file.exists(f)){next}
      A = kegg_list[[x[['A']]]]
      B = kegg_list[[x[['B']]]]
      
      df = data.frame(
        A = x[['A']],
        B = x[['B']],
        BMA = GOSemSim::clusterSim(A, B, semData=hsGO, measure="Wang", drop = "IEA", combine="BMA"),
        RCMAX =  GOSemSim::clusterSim(A, B, semData=hsGO, measure="Wang", drop = "IEA", combine="rcmax"),
        AVG =  GOSemSim::clusterSim(A, B, semData=hsGO, measure="Wang", drop = "IEA", combine="avg")
      )
      readr::write_tsv(df, file = f)
    }
  }
    
    
n = nrow(overlap_xtalk_df)
chunk_size = n %/% ncores
chunks = list()
for (i in 1:ncores) {
  start = (i - 1) * chunk_size + 1
  end = ifelse(i == ncores, n, i * chunk_size)
  # print( overlap_xtalk_df[start:end,])
  chunks[[i]] = overlap_xtalk_df[start:end,]
}

cl = parallel::makeCluster(ncores,outfile = 'GOsemantic_kegg.txt')
doSNOW::registerDoSNOW(cl)
# Apply the function to each chunk in parallel
boot = foreach::foreach(i = seq(ncores))
foreach::`%dopar%`(boot, compute_GOsemantic(chunks[[i]]))
parallel::stopCluster(cl)

## Gather all results into 1 single file 
for (ontology in c('BP','MF','CC')){
  df = apply(overlap_xtalk_df,1,function(x){
    read_tsv(file = paste0('data/KEGG/semantic/',ontology,'/',x[['A']],'-',x[['B']],'.tsv') ,col_types = cols())
  }) %>% data.table::rbindlist()
  
  write_tsv(df, file =  paste0('results/KEGG/GO-',ontology,'_SemanticSimilarity.tsv.gz'))
}
