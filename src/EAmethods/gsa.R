
run_gsa = function(se,
                     kegg_list,
                     seed_n){
  
  if(any(is.na(assay(se)))){
    ## GSA can't coope with NaN 
    se = subset(se,complete.cases(assay(se)))
  }
  
  res = 
    GSA::GSA(x=assay(se), 
             y=se[['GROUP']]+1, 
             nperms=2000, 
             genesets=kegg_list, 
             resp.type="Two class unpaired",
             genenames=names(se), 
             minsize=1,
             maxsize=1000,
             random.seed = seed_n)
  
  ps = cbind(res$pvalues.lo, res$pvalues.hi)
  ps = 2 * apply(ps, 1, min)
  scores = res$GSA.scores
  res = as.data.frame(cbind(scores, ps))
  colnames(res) = c("ES","Pvalue")

  return(res)
}
