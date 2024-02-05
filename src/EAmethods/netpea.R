
run_netpea = function(network_graph,
                      dea_genes,
                      kegg_df,
                      ncores,
                      seed_n){
  
  .enrichment_analysis = function(rwr_output_df){
    obs_col = 1
    perm_col = 2:2001
    
    tmp_enrichment = 
      data.frame(
        Ensembl = rownames(rwr_output_df),
        Scores = rwr_output_df[,obs_col]) %>%
      inner_join(kegg_df,by='Ensembl') %>%
      group_by(pathway_id) %>%

      summarise(Score = mean(Scores,na.rm=T)) %>%
      ungroup() %>%
      mutate(OBS=1)
    
    tmp_enrichment = 
      rbind(
        tmp_enrichment,
        lapply(perm_col,function(i){
          data.frame(
            Ensembl = rownames(rwr_output_df),
            Scores = rwr_output_df[,i]) %>%
            inner_join(kegg_df,by='Ensembl') %>%
            group_by(pathway_id) %>%
            summarise(Score = mean(Scores,na.rm=T)) %>%
            ungroup() %>%
            mutate(OBS=0)
        }) %>% data.table::rbindlist()
      )
    return(tmp_enrichment)
  }

  network_v = igraph::V(network_graph)$name
  
  restart_list = as.matrix(ifelse(network_v %in% dea_genes, 1, 0))
  rownames(restart_list) = network_v
  
  ### BUILD NULL MODEL
  permutations = 2000
  set.seed(seed_n)
  seeds_null = sample(10000,permutations,replace = F)
  nsize = length(dea_genes)
  to_sample = c(rep(1,nsize),rep(0,length(network_v)-nsize))
  
  restart_list = 
    cbind(
      restart_list,
        sapply(
          seq(permutations), function(i){
            set.seed(seeds_null[i])
            return(sample(to_sample,replace = F))
          })
      )
  
  set.seed(12345)
  rwr_output_df = 
    dnet::dRWR(network_graph,
               normalise = 'row',
               setSeeds = restart_list,
               restart = 0.5,
               verbose = FALSE,
               parallel = ncores>1,
               multicores = ncores
    )

  rownames(rwr_output_df) = network_v
  colnames(rwr_output_df) = c('OBS',paste0('P',seq(permutations)))
  tmp_enrichment_df = .enrichment_analysis(rwr_output_df)
  
  enrichment_df = tmp_enrichment_df %>%
    filter(OBS==0) %>%
    group_by(pathway_id) %>%
    summarise(avgScore = mean(Score,na.rm=T),
              sdScore = sd(Score,na.rm=T)) %>%
    left_join(tmp_enrichment_df %>%
                filter(OBS==1), by = 'pathway_id') %>%
    select(-OBS) %>%
    ungroup() %>%
    mutate(Z = (Score-avgScore)/sdScore) %>%
    mutate(Pvalue = pnorm(Z, mean = 0, sd = 1, lower.tail = FALSE))
  
  return(enrichment_df)
}
