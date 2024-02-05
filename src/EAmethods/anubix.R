
run_anubix = function(network,
                      network_matrix,
                      gs_df,
                      kegg_df,
                      ncores,
                      pfc,
                      nsamples,
                      network_type = 'weighted'){
  ANUBIX::anubix_constrained(
    network = network,
    links_matrix = network_matrix,
    genesets = gs_df,
    pathways = as.data.frame(kegg_df),
    cores=ncores,
    cutoff = pfc,
    sampling = nsamples,
    network_type = network_type,
    website = FALSE)
}