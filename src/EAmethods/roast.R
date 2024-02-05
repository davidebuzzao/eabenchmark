
run_roast = function(se,
                     kegg_list,
                     seed_n){
  
  group = factor(se$GROUP)
  design = model.matrix(formula('~group'))
  
  if(any(is.na(assay(se)))){
    ## ROAST can't coope with NaN
    se = subset(se,complete.cases(assay(se)))
  }
  
  set.seed(seed_n)
  limma::mroast(
    assay(se), 
    limma::ids2indices(kegg_list, rownames(se)), 
    design,
    nrot = 2000,
    adjust.method='none',
    sort='none'
  )
}