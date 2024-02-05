rm(list = ls())
setwd('~/eabenchmark/') 

library(gemma.R)
library(tidyverse)
library(data.table)
library(SummarizedExperiment)
library(parallel)

# Select case control studies
organism = 'hsa'
gemma.df = read_tsv(paste0('data/gemma/kegg_disease/',organism,'/summary.tsv'), col_types = 'ccdcdccl')
gemma_annotation.df = read_tsv(paste0('data/gemma/kegg_disease/',organism,'/annotation_summary.tsv'), col_types = 'cccc')

case_control.df = droplevels(
  gemma.df[gemma.df$ee.ShortName %in% 
             unique(
               as.character(
                 gemma_annotation.df[
                   grepl('disease',#phenotype
                         gemma_annotation.df$class.Name,
                         perl = T, 
                         ignore.case = T),
                   'ee.shortName'][[1]]))
           ,])
# length(unique(as.character(case_control.df$ee.ShortName))) #368 datasets

multiple_query_df = data.frame()
for (dataset in unique(case_control.df$ee.ShortName)){
  df = case_control.df[case_control.df$ee.ShortName==dataset,]
  tmp.df = df[1,]
  tmp.df$query = paste(df$query,collapse = '/')
  multiple_query_df = rbind(multiple_query_df,
                            tmp.df)
}

extract_col_an = function(d){
    x = multiple_query_df[multiple_query_df$ee.ShortName == 'd',]
    
    dataset = x[['ee.ShortName']]
    f.name = x[['ee.Name']]
    query = x[['query']]
    
    if (file.exists(paste0('data/gemma/',dataset,'.Rds'))){
      dat = readRDS(file = paste0('data/gemma/',dataset,'.Rds'))
    } else {
      dat <- tryCatch(
        {gemma.R::getDataset(dataset, 
                             filter = T,
                             type = "se")
        },
        error=function(cond) {
          return(NULL)})
      if (is.null(dat)){return()}
      
      saveRDS(dat,file = paste0('data/gemma/',dataset,'.Rds'))
      Sys.sleep(30)
    }
    
    data.frame(ee.ShortName = dataset,
               query = query,
               col_ano = paste(colnames(SummarizedExperiment::colData(dat)),collapse = '/'))
  }

ncores = 8
pb = utils::txtProgressBar(min=0, max=nrow(multiple_query_df), style = 3)
progress = function(n) utils::setTxtProgressBar(pb, n)
opts = list(progress = progress)
cl = parallel::makeCluster(ncores)
doSNOW::registerDoSNOW(cl)
boot = foreach::foreach(i = unique(multiple_query_df$ee.ShortName), .options.snow = opts)
gemma_df = foreach::`%dopar%`(boot, get_gemma_annotation(i)) %>% 
  data.table::rbindlist()
parallel::stopCluster(cl)


write_tsv(gemma_df, file = paste0('data/gemma/stats_',organism,'.tsv'))
