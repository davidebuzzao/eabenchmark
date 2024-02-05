#!/usr/bin/env Rscript

########################################################################
## Here you will find how to fetch human pathways info from KEGG-API
########################################################################
rm(list = ls())
setwd('~/eabenchmark/')
      
require(KEGGREST)
require(tidyverse)

## Choose the organism
organism = 'hsa' ##  22,243 entries (Release 101.0+/02-20, Feb 22)
# organism = 'mmu' ## 25,992 entries (Release 103.0+/08-09, Aug 22)
## displays the number of gene entries for the KEGG organism Homo sapiens
res = keggInfo(organism)
cat(res) 
## extract the latest list of KEGG pathways
databases = listDatabases()

query_tables = keggList('pathway', organism)
hsa_path = gsub('path:','', names(query_tables))

## extract genesID,genes-description 
## To get genes involved in a list of pathways
kegg.list = lapply(hsa_path,
                   function(pathway){
                     message(pathway)
                     res = keggGet(pathway)
                     if(!is.null(res[[1]]$GENE)){
                       data.frame(pathway_id = pathway,
                                  pathway_name = gsub(' - Homo sapiens \\(human\\)','',res[[1]]$NAME),
                                  # pathway_name = gsub(' - Mus musculus \\(house mouse\\)','',res[[1]]$NAME),
                                  entrezgene_id = res[[1]]$GENE[c(TRUE,FALSE)],
                                  hgnc_symbol = sapply(stringr::str_split(res[[1]]$GENE[c(FALSE,TRUE)],'; '), `[[`, 1),
                                  gene_description = res[[1]]$GENE[c(FALSE,TRUE)],
                                  class = stringr::str_split(res[[1]]$CLASS,'; ')[[1]][1],
                                  subclass = stringr::str_split(res[[1]]$CLASS,'; ')[[1]][2])
                     }
                   })

kegg.df = dplyr::bind_rows(kegg.list)
write_tsv(kegg.df,'data/KEGG/keggrest_hsa.tsv')