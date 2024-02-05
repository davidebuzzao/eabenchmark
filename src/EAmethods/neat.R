
run_neat = function(network,
                    network_v,
                    gs_list,
                    kegg_list
                    ){
  neat::neat(
    alist=gs_list, 
    blist=kegg_list,
    network = network,
    nodes = network_v,
    nettype = 'undirected',
    mtc.type = 'none')
}