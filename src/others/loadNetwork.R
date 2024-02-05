
library(igraph)

load_network = function(network.name,cutoff=0.8){
  if (network.name == 'funcoup'){
    network = read_tsv('data/network/FC5.0_H.sapiens_compact.gz', col_types = cols())
    network = network[,c('2:Gene1', '3:Gene2','#0:PFC')]; colnames(network) = c('gene1','gene2','PFC')
    network = network[network$PFC>=cutoff, c('gene1', 'gene2')]
    
    graph = igraph::graph_from_data_frame(network, directed = FALSE)
    
  } else if (network.name == 'string'){
    ## Load STRING, simplify the graph and filter out low quality links
    version = '11.5'
    network = read_csv(paste0('data/network/9606.protein.links.v',version,'.txt.gz'), col_types = cols())
    ## 
    network$combined_score = network$combined_score/1000
    network.df = network[network$combined_score>=cutoff,]
    graph = igraph::graph_from_data_frame(network.df, directed = FALSE)
  } else if (network.name == 'humanet'){
    network = dplyr::tibble(read.table('data/network/HumanNet-XC.tsv.gz', sep = '', header = FALSE, col.names = c('EntrezGeneID1','EntrezGeneID2','LLS')))
    network = network[order(network$LLS, decreasing=T),]
    if(cutoff=='10p'){
      network = network[1:round(nrow(network)/10),]
    } else if (cutoff=='20p'){
      network = network[1:round(nrow(network)/5),] #225099
    }

    graph = igraph::graph_from_data_frame(network, directed = FALSE)
    
  } 
  
  ## Remove duplicated links
  if (!igraph::is_simple(graph)){
    graph = igraph::simplify(graph,remove.multiple = TRUE, remove.loops = TRUE)
  }
  
  return(graph)
}
