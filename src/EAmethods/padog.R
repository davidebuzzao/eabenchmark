
run_padog = function(se,
                   kegg_list,
                   ncores,
                   seed_n){
  PADOG::padog(
    esetm = as.matrix(assay(se)),
    group = ifelse(se$GROUP==0,'c','d'),
    paired = FALSE,
    block = NULL,
    gslist = kegg_list,
    Nmin = 1,
    NI = 2000, 
    verbose = F,
    parallel = TRUE,
    ncr = ncores,
    dseed = seed_n
  )

}
