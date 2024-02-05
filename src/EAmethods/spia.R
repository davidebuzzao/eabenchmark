
run_spia = function(se,
                    GeneID,
                    beta_v,
                    seed_n){
  dea = as.data.frame(SummarizedExperiment::rowData(se)) 

  AllGenes = rownames(dea)
  DEGenes_logFC = as.double(dea[GeneID,'logFC'])
  names(DEGenes_logFC) = as.character(GeneID)
  
  set.seed(seed_n)
  SPIA::spia(de=DEGenes_logFC, 
             all=AllGenes, 
             beta = beta_v,
             nB = 2000,
             organism="hsa",
             data.dir="data/spia/",
             verbose = T)
}
