rm(list = ls())
setwd('~/eabenchmark/') 

library(tidyverse)
library(SummarizedExperiment)

## FUNCOUP mapping file
uniprot_mapping_df = 
  read_tsv('data/mapping/UP000005640_9606.idmapping.gz', 
           col_names = c('Uniprot_id','db','db_id'), 
           col_types = 'cfc') %>%
  mutate(db_id=gsub('hsa:','',db_id)) 

human_geneid = 
  uniprot_mapping_df %>%
  filter(db=='GeneID') %>%
  select(Uniprot_id,db_id) %>% 
  dplyr::rename(GeneID=db_id) %>%
  pull(GeneID) %>% unique()


cc_design_step3.df = read_tsv('data/TPbenchmark/cc_design_step3.tsv') 

cc_design_step3.df = cc_design_step3.df %>% 
  filter(control.Samples>2 & 
           disease.Samples >2 & 
           !grepl('\\.[0-9]',ee.ShortName,perl = T) & 
           query!='Primary immunodeficiency' & 
           ! ee.Accession %in% c('GSE24283','GSE5967','GSE26637', 'GSE4302', 'GSE11543') & 
           ## batch effect + bad normalization
           ! ee.Accession %in% c('GSE48301','GSE44000'))

cc_design_step3.df = lapply(unique(cc_design_step3.df$ee.ShortName), function(x){
  df = cc_design_step3.df %>% filter(ee.ShortName == x)
  control = unique(df$control)
  df$control = paste(control,collapse = ',')
  disease = unique(df$disease)
  df$disease = paste(disease,collapse = ',')
  return(unique(df))
}) %>% data.table::rbindlist()

dataset_df = cc_design_step3.df %>%
  dplyr::rename(TargetPathway = query,
                GEO = ee.ShortName, 
                Title = ee.Name,
                DataType = data.type,
                BatchEffect = ee.batchEffect) %>%
  select(GEO,TargetPathway,Title,DataType,BatchEffect)

logtransform_data = function(edata, ano_col){
  ## GEO2R scripts; source biostar: https://www.biostars.org/p/99879/
  ex = edata
  qx = as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC = (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) { 
    ex[which(ex <= 0)] = NaN
    return(log2(ex[,ano_col]))
  } else {return(ex[,ano_col])}
}

alpha = c(0.1, 0.2)
min_coverage = 1e4
max_n = 500
min_n = 15
count = 1; tot = nrow(dataset_df)

out_df = apply(dataset_df, 1, function(x){
  message(paste0(x[['GEO']],' ', count,'/',tot)); count <<- count +1
  if (!file.exists(paste0('data/TPbenchmark/',x[['DataType']],'/',x[['GEO']],'.Rds'))){return()}
  se = readRDS(paste0('data/TPbenchmark/',x[['DataType']],'/',x[['GEO']],'.Rds'))
  if (x[['DataType']]=='rseq'){
    se = se[rownames(se) %in% human_geneid]
    saveRDS(se,paste0('data/TPbenchmark/',x[['DataType']],'/',x[['GEO']],'.Rds'))
  }
  tmp_dat = assay(se)
  ano = colData(se)
  dat = as.matrix(logtransform_data(tmp_dat,rownames(ano)))
  if(!all(tmp_dat == dat, na.rm = T)){
    message(paste0(dataset, ' --> log-transform'))
  }
  if(!all(tmp_dat == dat, na.rm = T)){
    message(paste0(dataset, ' --> log-transform'))
    f = '~group'
    group = factor(ano$GROUP)
    f = formula(f)
    design = model.matrix(stats::formula(f))

    fit = limma::lmFit(dat, design)

    if (data.type == 'ma'){
      fit = limma::eBayes(fit)
    } else if (data.type == 'rseq'){
      fit = limma::eBayes(fit,trend = T)
    }

    aT1 = limma::topTable(fit, number=nrow(dat), coef="group1",
                          sort.by="none",adjust.method="BH")
    de = aT1[, c("logFC", "P.Value", "t", "adj.P.Val")]
    as.data.frame(de) %>% arrange(adj.P.Val) %>% head

    # create new SE
    newDat = SummarizedExperiment(assays=dat,
                                  rowData = de,
                                  colData=colData(se),
                                  metadata=metadata(se))

    saveRDS(newDat,file = paste0('data/TPbenchmark/',data.type,'/',dataset,'.Rds'))

    se = readRDS(paste0('data/TPbenchmark/',data.type,'/',dataset,'.Rds'))
    dat = assay(se)
  }
  
  genome_coverage = nrow(dat)
  dea_df = as.data.frame(rowData(se)) 
  
  alpha = 0.1
  n_deg = nrow(
    dea_df %>%
      filter(adj.P.Val < alpha)
  )
  
  if (n_deg>max_n){
    n_deg = max_n
  } else if (n_deg<min_n){
    alpha = 0.2
    n_deg = nrow(
      dea_df %>%
        filter(adj.P.Val < alpha)
    )
    if (n_deg>max_n){
      n_deg = max_n
    } else if (n_deg<min_n){
      n_deg = NA
    }
  }
  
  out_df = 
    data.frame(
      GEO = x[['GEO']],
      Title = x[['Title']],
      TargetPathway = x[['TargetPathway']],
      DataType = x[['DataType']],
      BatchEffect = x[['BatchEffect']],
      GenomeCoverage = genome_coverage,
      NrOfDEG = n_deg,
      rankBy = 't-statistics',
      alpha = alpha,
      beta = 0,
      MinMax = paste(min_n,max_n,sep='/')
    )
  
  return(out_df)
}) %>% data.table::rbindlist()

write_tsv(out_df,file=paste0('data/TPbenchmark/cc_deg-alpha10OR20-Min',min_n,'Max',max_n,'-beta0.tsv'))
