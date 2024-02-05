#!/usr/bin/env Rscript

run_ease = function(A,PT=18000,gs){                          
  # https://david.ncifcrf.gov/helps/functional_annotation.html
  # LH: list hits
  # LT: list total
  # PH: population hits
  # PT: population total
  # User genes | genome
  # ------------------
  # LH-1    | PH-LH+1 | PH
  # ------------------
  # LT-LH   | PT-LT-(PH-LH) | PT-PH
  # ------------------
  
  A = as.character(A)
  LT = length(A)
  PT = as.numeric(PT) # experimental
  
  gs_names = names(gs)
  enrichment_df = lapply(gs_names, function(gs_n){
    B = as.character(gs[[gs_n]])
    PH = length(B)
    LH = length(intersect(A,B))
    
    modLH = ifelse(LH==0,1,LH)
    ## Implement EASEscore from DAVID 
    contingency_table = rbind(c(modLH-1, PH-LH+1),
                              c(LT-LH,PT-LT-(PH-LH))) 
    
    pval=fisher.test(contingency_table, alternative = 'greater')[['p.value']]
    
    return(
      data.frame(
        pathway_id = gs_n,
        Overlap = LH,
        Pvalue = pval,
        GeneID = ifelse(LH>0,paste0(intersect(A,B), collapse = ','),NA))
    )   
  }) %>% data.table::rbindlist()
  return(enrichment_df)
}
