
run_cepa = function(se,
                    GeneID,
                    kegg_catalogue,
                    method = 'ora',
                    seed_n){
  ## use all 6 pathways centralities, get the lowest pval
  library(CePa)
  if(method =='ora'){
    set.seed(seed_n)
    cepa.all(dif = GeneID, 
             bk = rownames(se), 
             pc = kegg_catalogue, 
             iter = 2000)
  } else if (method =='gsa'){
    set.seed(seed_n)
    cepa.all(mat = as.matrix(assay(se)), 
             label = sampleLabel(se$GROUP, treatment = 1, control = 0), 
             pc = kegg_catalogue, 
             iter = 2000)
  }
}