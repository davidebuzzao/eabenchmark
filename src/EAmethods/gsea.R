
run_gsea = function(se,
                   kegg_df,
                   seed_n){
  
  dea_df = as.data.frame(SummarizedExperiment::rowData(se)) 
  dea_df = dea_df[order(dea_df$t, decreasing = T),]
  dea_ranks = dea_df$t
  names(dea_ranks) = rownames(dea_df)
  dea_ranks = dea_ranks[!is.na(dea_ranks)]
  
  set.seed(seed_n)
  clusterProfiler::GSEA(
    dea_ranks,
    nPerm=2000,
    exponent = 1,
    minGSSize = 1,
    maxGSSize = 1000,
    eps = 1e-10,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    TERM2GENE = as.data.frame(kegg_df),
    TERM2NAME = NA,
    verbose = FALSE,
    seed = seed_n,
    by = "DOSE",
  )
}

