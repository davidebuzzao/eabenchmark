
library(tidyverse)

load_mappingFile = function(network.name){
  
  if (network.name == 'funcoup'){
    uniprot_mapping_df = 
      read_tsv('data/mapping/UP000005640_9606.idmapping.gz', 
               col_names = c('Uniprot_id','db','db_id'), 
               col_types = 'cfc') %>%
      mutate(db_id=gsub('hsa:','',db_id)) %>%
      filter(db %in% c('Ensembl','GeneID'))
    
    return(uniprot_mapping_df %>%
            filter(db=='GeneID') %>%
            dplyr::rename(GeneID=db_id) %>%
             ## ENSEMBL
            left_join(uniprot_mapping_df %>%
                  filter(db=='Ensembl') %>%
                  dplyr::rename(Ensembl=db_id), by='Uniprot_id', multiple='all') %>%
             select(GeneID,Ensembl) %>% 
             unique() 
    )
    
  } else if (network.name == 'string'){
    string_mapping_df = 
      read_tsv('data/mapping/9606.protein.aliases.v11.5.txt.gz', 
               col_names = c('string_protein_id','alias','source'), 
               col_types = 'ccf', comment = '#') %>%
      filter(grepl('Gene_ID|GeneID',source, perl = T)) %>%
      mutate(string_protein_id = gsub('9606\\.','',string_protein_id)) %>%
      dplyr::rename(Ensembl=string_protein_id,
                    GeneID=alias)
    
    return(
      string_mapping_df %>% 
        select(GeneID, Ensembl) %>%
        unique()
    )
  }
}