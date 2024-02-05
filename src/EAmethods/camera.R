
run_camera = function(se,
                      kegg_list,
                      var_corr = F){

  group = factor(se$GROUP)
  design = model.matrix(formula('~group'))
  
  if(any(is.na(assay(se)))){
    ## Camera can't coope with NaN
    se = subset(se,complete.cases(assay(se)))
  }
  
  ## trend.var --> TRUE if rseq, FALSE if ma
  if (!var_corr){
    #camera will rank biologically interpretable sets more highly. 
    # This gives a useful compromise between strict error rate control and interpretable gene set rankings.
    limma::camera(
      assay(se), 
      limma::ids2indices(kegg_list, rownames(se)), 
      design, 
      contrast = ncol(design),
      inter.gene.cor = 0.01,
      trend.var = se@metadata$dataType == 'ma',
      sort=FALSE
    )
  } else {
    ## an inter-gene correlation will be estimated for each tested set.
    limma::camera(
      assay(se), 
      limma::ids2indices(kegg_list, rownames(se)), 
      design, 
      contrast = ncol(design),
      inter.gene.cor = NA, 
      trend.var = se@metadata$dataType == 'ma',
      sort=FALSE
    )
  } 
}