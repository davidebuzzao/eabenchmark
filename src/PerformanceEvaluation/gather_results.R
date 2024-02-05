rm(list = ls())
setwd('~/eabenchmark/') 

library(tidyverse)

networks = c('funcoup','string')
ora = c('ease','fisher')
anubix = paste('anubix',networks,sep = '-')
neat = paste('neat',networks,sep = '-')
binox = paste('binox',networks,sep = '-')
netpea = paste('netpea',networks,sep = '-')
camera = c('camera','camera*')
methods = c(anubix,binox,neat,netpea,camera,ora,'roast','fgsea','gsa','gsea','padog','gsva','cepaORA','spia')

## KEGG
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

## The expression data
dataset_df = read_tsv(file='data/TPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv', 
                      col_types = 'fcfffiddiciccccc')

## The disease pathway network
top_cutoff = c(10,20,40,'ALL')
disease_network = lapply(top_cutoff, function(top){
  tmp_df = read_tsv(file = paste0('results/DiseasePathwayNetwork/DiseasePathwayNetwork-FDR<5e-2_top',top,'GOaverageBMA.tsv.gz'), col_types = cols())
  tmp_df$GO_BMA = round(tmp_df$GO_BMA,3)
  tmp_df$TOP = top
  return(tmp_df)
}) %>% data.table::rbindlist()


parse_results = function(method, benchmark='TPbenchmark',fill_missing_pathway=T){

  if (benchmark == 'TPbenchmark'){
    res = read_tsv(file=paste0('results/TPbenchmark/',method,'.tsv.gz'), col_types = cols()) %>%
      dplyr::filter(gsub('\\.[0-9]*','',GEO) %in% unique(dataset_df$GEO) &
                      pathway_id %in% unique(kegg_df$pathway_id))
  } else if (benchmark == 'FPbenchmark'){
    res = read_tsv(file=paste0('results/FPbenchmark/',method,'.tsv.gz'), col_types = cols()) %>%
      dplyr::filter(gsub('\\.[0-9]*','',GEO) %in% unique(dataset_df$GEO) &
                      pathway_id %in% unique(kegg_df$pathway_id))
  } 

  if (method %in% c(binox,neat)){
    ## Treat depletions a poor enrichments
    res = res %>%
      mutate(Pvalue = ifelse(nab < expected_nab, 1.0 - as.double(Pvalue), Pvalue))
  } else if (method %in% 'cepaORA'){
    ## Extract min pval
    res = res %>%
      group_by(GEO,pathway_id) %>%
      summarise(Pvalue = min(equal.weight,in.degree,
                             out.degree,betweenness,
                             in.reach,out.reach)) %>%
      ungroup()
  } 
  
  res = res %>%
    select(GEO,pathway_id,Pvalue)
  
  if(fill_missing_pathway & benchmark == 'TPbenchmark'){
    kegg_pathways = unique(kegg_df$pathway_id)
    res = lapply(unique(res$GEO), function(geo){
      tmp_res = res %>%
        dplyr::filter(GEO == geo) %>% droplevels()
      missing_pathways = setdiff(kegg_pathways,tmp_res$pathway_id)
      if (length(missing_pathways)==0) return(tmp_res)
      return(rbind(tmp_res,
                   data.frame(GEO = geo,
                              pathway_id = missing_pathways,
                              Pvalue = 1)))
    }) %>% data.table::rbindlist()
    
  } else if(fill_missing_pathway & benchmark == 'FPbenchmark'){
    kegg_pathways = unique(kegg_df$pathway_id)
    fill_pathway = 
      function(geo){
        library(tidyverse)
        tmp_res = res %>%
          dplyr::filter(GEO == geo) %>% droplevels()
        missing_pathways = setdiff(kegg_pathways,tmp_res$pathway_id)
        if (length(missing_pathways)==0) return(tmp_res)
        return(rbind(tmp_res,
                     data.frame(GEO = geo,
                                pathway_id = missing_pathways,
                                Pvalue = 1)))
      }
    
    pb = utils::txtProgressBar(min=0, max=length(unique(res$GEO)), style = 3)
    progress = function(n) utils::setTxtProgressBar(pb, n)
    opts = list(progress = progress)
    cl = parallel::makeCluster(8)
    doSNOW::registerDoSNOW(cl)
    boot = foreach::foreach(i = unique(res$GEO), .options.snow = opts)
    res = foreach::`%dopar%`(boot, fill_pathway(i)) %>% 
      data.table::rbindlist()
    parallel::stopCluster(cl)
  } 
  return(res)
}


########################
## TP benchmark
positive_benchmark = lapply(methods, function(m){
  message(m)
  tmp_df = parse_results(m, benchmark='TPbenchmark', fill_missing_pathway=T)
  nd = setdiff(unique(dataset_df$GEO),unique(tmp_df$GEO))
  if (length(nd)>0){
    print(paste0(m,' has 0 predictions on ', length(nd), ' datasets.'))
    kegg_pathways = unique(kegg_df$pathway_id)
    tmp_df = bind_rows(tmp_df,
                       lapply(nd, function(d){
                         data.frame(GEO=d,
                                    pathway_id=kegg_pathways,
                                    Pvalue = 1)
                       }) %>% data.table::rbindlist()
    )
  }
  tmp_df$Method = m
  
  enrichments_df = lapply(unique(tmp_df$GEO), function(geo){
    target = as.character(
      dataset_df %>%
        dplyr::filter(GEO==geo) %>% droplevels() %>%
        pull(TargetPathway))
    
    df = lapply(top_cutoff, function(t){
      tmp2_df = tmp_df %>%
        dplyr::filter(GEO==geo) %>% droplevels()
      
      tmp_disease_network = disease_network %>% 
        dplyr::filter(TOP==t & A_name==target)
      
      tmp2_df$label = ifelse(tmp2_df$pathway_id %in% unique(c(tmp_disease_network$A,tmp_disease_network$B)),1,0)
      tmp2_df$TOP = t
      return(tmp2_df)
      }) %>% data.table::rbindlist()
    return(df)
  }) %>% data.table::rbindlist()
  
  return(enrichments_df)
}) %>%  data.table::rbindlist()

positive_benchmark %>%
  group_by(Method) %>%
  summarise(n())

write_tsv(positive_benchmark, file = 'results/stats/TPbenchmark.tsv.gz')
# write_tsv(positive_benchmark %>% filter(label==1), file = 'results/stats/TPbenchmark_only1.tsv.gz')

########################
## FP benchmark - unbalanced
negative_dataset_df = read_tsv(file = 'data/FPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv.gz',
                               col_types = 'fccffiddiciccccci') 

negative_benchmark = lapply(methods, function(m){
  message(paste0('\n',m))
  tmp_df = parse_results(m, benchmark='FPbenchmark', fill_missing_pathway=T)
  nd = setdiff(unique(negative_dataset_df$GEO),unique(tmp_df$GEO))
  if (length(nd)>0){
    print(paste0(m,' has 0 predictions on ', length(nd), ' datasets.'))
    kegg_pathways = unique(kegg_df$pathway_id)
    tmp_df = bind_rows(tmp_df,
                       lapply(nd, function(d){
                         data.frame(GEO=d,
                                    pathway_id=kegg_pathways,
                                    Pvalue = 1)
                       }) %>% data.table::rbindlist()
    )
  }
  tmp_df$Method = m
  return(tmp_df)
  
}) %>% data.table::rbindlist()

negative_benchmark %>%
  group_by(Method) %>%
  summarise(n()) %>% View()

write_tsv(negative_benchmark, file = 'results/stats/FPbenchmark.tsv.gz')

########################
## FP benchmark - balanced

## Only few pathways per dataset is negative
negative_benchmark = lapply(methods, function(m){
    message(m)
    ##### subset data
    tmp_dataset_df = dataset_df %>%
      dplyr::select(GEO,TargetPathway)
    
    tmp_df = parse_results(m, benchmark='FPbenchmark', fill_missing_pathway=T)
    nd = setdiff(unique(negative_dataset_df$GEO),unique(tmp_df$GEO))
    if (length(nd)>0){
      print(paste0(m,' has 0 predictions on ', length(nd), ' datasets.'))
      kegg_pathways = unique(kegg_df$pathway_id)
      tmp_df = bind_rows(tmp_df,
        lapply(nd, function(d){
        data.frame(GEO=d,
                   pathway_id=kegg_pathways,
                   Pvalue = 1)
          }) %>% data.table::rbindlist()
      )
    }
    tmp_df$Method = m
    
    extract_tp = 
      function(geo){
        library(tidyverse)
        top_cutoff = c(10,20,40,'ALL')
        disease_network = lapply(top_cutoff, function(top){
          tmp_df = read_tsv(file = paste0('results/DiseasePathwayNetwork/DiseasePathwayNetwork-FDR<5e-2_top',top,'GOaverageBMA.tsv.gz'), col_types = cols())
          tmp_df$TOP = top
          return(tmp_df)
        }) %>% data.table::rbindlist()
        
        tmp2_df = tmp_df %>%
          dplyr::filter(GEO==geo) %>% droplevels()
        target = as.character(
          tmp_dataset_df %>%
            dplyr::filter(GEO==gsub('\\.[0-9]*','',geo)) %>% droplevels() %>%
            pull(TargetPathway))
        
        out_df = lapply(top_cutoff, function(t){
          tmp_disease_network = disease_network %>% 
          filter(TOP==t) %>% droplevels()
          
          p = tmp_disease_network %>% 
           dplyr::filter(A_name==target) %>% droplevels()
          
          tmp2_df$label = ifelse(tmp2_df$pathway_id %in% unique(c(p$A,p$B)),1,0)
          tmp2_df$TOP = t
          return(tmp2_df)
          
        }) %>% data.table::rbindlist()
        return(out_df)
        }
    
    pb = utils::txtProgressBar(min=0, max=length(unique(tmp_df$GEO)), style = 3)
    progress = function(n) utils::setTxtProgressBar(pb, n)
    opts = list(progress = progress)
    cl = parallel::makeCluster(10)
    doSNOW::registerDoSNOW(cl)
    boot = foreach::foreach(i = unique(tmp_df$GEO), .options.snow = opts)
    out_df = foreach::`%dopar%`(boot, extract_tp(i)) %>% 
      data.table::rbindlist()
    parallel::stopCluster(cl)
    
    return(out_df)
    }) %>% data.table::rbindlist()

negative_benchmark %>%
  select(Method,GEO) %>% unique() %>% 
  group_by(Method) %>%
  summarise(n()) 

write_tsv(negative_benchmark, file = 'results/stats/FPbenchmark_balanced.tsv.gz')

#### Count min,max,avg number of enrichments
positive_benchmark_count = lapply(methods, function(m){
  message(m)
  tmp_df = parse_results(m, benchmark='positive', fill_missing_pathway=F)
  tmp_df = tmp_df %>% 
    group_by(GEO) %>%
    summarise(tot_pathways=n()) 
  
  nd = setdiff(unique(dataset_df$GEO),unique(tmp_df$GEO))
  if (length(nd)>0){
    print(paste0(m,' has 0 predictions on ', length(nd), ' datasets.'))
    kegg_pathways = unique(kegg_df$pathway_id)
    tmp_df = bind_rows(tmp_df,
                       lapply(nd, function(d){
                         data.frame(GEO=d,
                                    tot_pathways = 0)
                       }) %>% data.table::rbindlist()
    )
  }
  
  return(
    data.frame(
      Method = m,
      val = c(min(tmp_df$tot_pathways),max(tmp_df$tot_pathways), median(tmp_df$tot_pathways)),
      type = c('min','max','median'))
  )
}) %>%  data.table::rbindlist()

write_tsv(positive_benchmark_count, file = 'results/stats/TPbenchmark_counts.tsv')

