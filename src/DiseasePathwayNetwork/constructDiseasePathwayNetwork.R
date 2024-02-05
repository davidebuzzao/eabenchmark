### Both overlap and network connectivity have been tested bidirectionally with degree-aware resampling.
### Here we combine (Fisher method) pvalues and extract FDR independently for overlap and network separation.
## Later, we take the intersection of significant results as disease pathway network
rm(list = ls())
setwd('~/eabenchmark/') 

library(tidyverse)
library(igraph)
library(metap)

#######################################
## Disease Pathway Network
#######################################

## Inter-pathway overlap
overlap_df = 
  read_tsv('results/KEGG/InterPahtwayOverlap_degree-aware1kSampling_KEGGinHumanNet20p.tsv.gz', 
           col_select = c('A','B','direction','overlap','sample_pval'), 
           col_types = 'ccfdd') %>%
  dplyr::rename(overlap_pval = sample_pval) %>%
  # Swap and sort columns
  rowwise() %>%
  mutate(A = sort(c(A, B))[1],
         B = sort(c(A, B))[2]) %>%
  ungroup()

## Inter-pathway connectivity
xtalk_df = read_tsv('results/KEGG/InterPahtwayConnectivity_degree-aware1kSampling_KEGGinHumanNet20p.tsv.gz', 
                    col_select = c('A','B','direction','Xtalk','sample_pval'), 
                    col_types = 'ccfdd') %>%
  dplyr::rename(xtalk_pval = sample_pval) %>%
  # Swap and sort columns
  rowwise() %>%
  mutate(A = sort(c(A, B))[1],
         B = sort(c(A, B))[2]) %>%
  ungroup()

overlap_xtalk_df = overlap_df %>% 
  inner_join(xtalk_df, by = c('A','B','direction'), multiple = 'all') %>% 
  group_by(A,B) %>% 
  summarise(overlap_sumlog = metap::sumlog(overlap_pval)[['p']],
            xtalk_sumlog = metap::sumlog(xtalk_pval)[['p']]) %>%
  ungroup() %>%
  mutate(overlap_FDR = p.adjust(overlap_sumlog, method = 'BH'),
         xtalk_FDR = p.adjust(xtalk_sumlog, method = 'BH')) %>%
  filter(overlap_FDR<0.05 & xtalk_FDR<0.05)

avg_overlap_xtalk_df = overlap_xtalk_df %>% 
  group_by(A) %>% 
  count()

##################################################
########## GO SEMANTIC SIMILARITY
ontologies = c('BP','MF','CC')
SemanticSimilarity_df = lapply(ontologies, function(ont){
  sem_sim_df = 
    read_tsv(file = paste0('results/KEGG/GO-',ont,'_SemanticSimilarity.tsv.gz'), skip = 1,
             col_names =  c('A','B',paste0(ont,'_BMA'),paste0(ont,'_RCMAX'),paste0(ont,'_AVG')),
             col_types = 'ccddd') 
  return(sem_sim_df)
})

semantic_df = SemanticSimilarity_df[[1]] %>% 
  inner_join(SemanticSimilarity_df[[2]], by = c('A','B')) %>% 
  inner_join(SemanticSimilarity_df[[3]], by = c('A','B')) %>% 
  rowwise() %>% 
  mutate(GO_BMA = mean(c(BP_BMA,MF_BMA,CC_BMA)),
         GO_RCMAX = mean(c(BP_RCMAX,MF_RCMAX,CC_RCMAX)),
         GO_AVG = mean(c(BP_AVG,MF_AVG,CC_AVG))) %>% 
  ungroup() %>% 
  select(A,B,GO_BMA,GO_RCMAX,GO_AVG)
write_tsv(semantic_df, file = 'results/KEGG/GOaverageBMA_SemanticSimilarity.tsv.gz')


##################################################
########## EXTRACT DISEASE NETWORK
semantic_df = read_tsv( file = 'results/KEGG/GOaverageBMA_SemanticSimilarity.tsv.gz')

pathway_network_df = overlap_xtalk_df %>%
  left_join(semantic_df, by = c('A','B')) %>% 
  replace(is.na(.),0)

dataset_df = read_tsv(file='data/TPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv', 
                      col_types = 'fcfffiddiciccccc') %>% 
  select(GEO,TargetPathway,) %>%
  droplevels()

## Load KEGG
kegg_df = read_tsv(file = 'data/KEGG/keggrest_hsa.tsv', col_types = 'ffcccff') 
all_pathway = kegg_df %>% 
  filter(pathway_id %in%
           as.character(
             kegg_df %>%
               select(pathway_id,entrezgene_id) %>%
               unique() %>%
               group_by(pathway_id) %>%
               summarise(n=n()) %>%
               filter(n>14) %>%
               droplevels() %>%
               pull(pathway_id))) %>%
  droplevels() %>%
  pull(pathway_id) %>% as.character() %>% unique()

kegg_human_df = kegg_df %>% 
  filter(class=='Human Diseases') %>% 
  select(pathway_id,pathway_name) %>% 
  unique()

pathway_disease = kegg_df %>%
  filter(pathway_name %in% dataset_df$TargetPathway) %>% 
  pull(pathway_id) %>% as.character() %>% unique()

pathway_network_df = pathway_network_df %>% 
  filter(A %in% all_pathway & B %in% all_pathway)

pathway_graph = igraph::graph_from_data_frame(pathway_network_df[,c('A','B','GO_BMA')], directed = FALSE)

pathway_graph
if (!is_simple(pathway_graph)){
  pathway_graph = simplify(pathway_graph, remove.multiple = TRUE, remove.loops = TRUE)
  # pathway_graph = simplify(pathway_graph, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = 'max')
}

message(paste('Graph loaded\n',
              '# Nodes: ',vcount(pathway_graph),'\n',
              '# Edges: ',ecount(pathway_graph),sep = ''))

## Extract TOP scoring links
d = degree(pathway_graph)

##################### 
# Extract TOP
## Load KEGG
kegg_df = read_tsv(file = 'data/KEGG/keggrest_hsa.tsv', col_types = 'ffcccff') %>%
  select(pathway_id,pathway_name,class,subclass) %>%
  unique() %>%
  droplevels()

top_cutoff = c(10,20,40,-1)
for (top in top_cutoff){
  message(top)
  pathwayome_topX_df = lapply(as.character(kegg_human_df$pathway_id),function(path){
    g = E(pathway_graph)[.from(path)]
    n = length(g)
    if (n<top){
      top = n
      print(paste("Switching to ",top))
    }
    if (top==0){return()}
    else if (top==-1){
      ## TAKE ALL
      tmp_df = data.frame( ids = as_ids(g)) 
      tmp_df = apply(tmp_df,1,function(x){
        p = unlist(stringr::str_split(x,'\\|'))
        data.frame(A=path,
                   B=p[!p %in% path])
      }) %>% data.table::rbindlist()
      tmp_df$GO_BMA = apply(tmp_df,1,function(x){
        E(pathway_graph,path=c(x[['A']],x[['B']]))$GO_BMA
      })
      return(tmp_df)
    } else {
      min_val = sort(g$GO_BMA,decreasing=T)[1:top][top]
      g = g[!g$GO_BMA<min_val]
      tmp_df = data.frame( ids = as_ids(g)) 
      
      tmp_df = apply(tmp_df,1,function(x){
        p = unlist(stringr::str_split(x,'\\|'))
        data.frame(A=path,
                   B=p[!p %in% path])
      }) %>% data.table::rbindlist()
      tmp_df$GO_BMA = apply(tmp_df,1,function(x){
        E(pathway_graph,path=c(x[['A']],x[['B']]))$GO_BMA
      })
      
      return(tmp_df)
    }
    
  }) %>% data.table::rbindlist()
  
  pathwayome_topX_df = pathwayome_topX_df %>% 
    left_join(kegg_df %>% 
                select(pathway_id,pathway_name) %>% 
                dplyr::rename(A=pathway_id),
              by = 'A') %>% 
    dplyr::rename(A_name=pathway_name) %>% 
    left_join(kegg_df %>% 
                select(pathway_id,pathway_name) %>% 
                dplyr::rename(B=pathway_id),
              by = 'B') %>% 
    dplyr::rename(B_name=pathway_name)
  print(pathwayome_topX_df)
  write_tsv(pathwayome_topX_df, file = paste0('results/DiseasePathwayNetwork/DiseasePathwayNetwork-FDR<5e-2_top',ifelse(top!=-1,top,'ALL'),'GOaverageBMA.tsv.gz'))
}

## Output table figure
## Load KEGG
kegg_df = read_tsv(file = 'data/KEGG/keggrest_hsa.tsv', col_types = 'ffcccff')
kegg_df = kegg_df %>%
  filter(pathway_id %in%
           as.character(
             kegg_df %>%
               select(pathway_id,entrezgene_id) %>%
               unique() %>%
               group_by(pathway_id) %>%
               summarise(n=n()) %>%
               filter(n>14) %>%
               droplevels() %>%
               pull(pathway_id))) %>%
  droplevels()

dataset_df = read_tsv(file='data/TPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv', 
                      col_types = 'fcfffiddiciccccc') %>% 
  select(GEO,TargetPathway,) %>%
  droplevels()

# Xs = c(1,5,10,25,50,100)
top_cutoff = c(10,20,40,'ALL')
pathway_network_df = lapply(top_cutoff, function(top){
  tmp_df = read_tsv(file = paste0('results/DiseasePathwayNetwork/DiseasePathwayNetwork-FDR<5e-2_top',ifelse(top!=-1,top,'ALL'),'GOaverageBMA.tsv.gz'), col_types = cols())
  
  tmp_df$TOP = top
  return(tmp_df)
}) %>% data.table::rbindlist()

pathway_network_df %>% 
  group_by(A,A_name,TOP) %>%
  dplyr::count(A) %>% 
  rowwise() %>%
  mutate(n=n+1) %>%
  ungroup() %>%
  pivot_wider(names_from = TOP, values_from = n) %>%
  replace(is.na(.), 1) %>%
  dplyr::rename(pathway_name = A_name,
                pathway_id = A) %>%
  write_tsv('notebooks/tables/DiseasePathwayNetwork_TOP10-40.tsv')

d = pathway_network_df %>% 
  group_by(A,A_name,TOP) %>%
  dplyr::count(A) %>%
  ungroup() %>%
  pull(n) %>% 
  mean()

pathway_network_df %>% 
  filter(A_name %in% dataset_df$TargetPathway) %>%
  group_by(A,A_name,TOP) %>%
  dplyr::count(A) %>% 
  rowwise() %>%
  mutate(n=n+1) %>%
  ungroup() %>%
  pivot_wider(names_from = TOP, values_from = n) %>%
  replace(is.na(.), 1) %>%
  dplyr::rename(pathway_name = A_name,
                pathway_id = A) %>%
  write_tsv('notebooks/tables/TargetDiseasePathwayNetwork_TOP10-ALL.tsv')


################# How frequently pathways are tested?
## Logic: For each dataset, list the tested pathways [target + top20]. Later, extract frequencies.
dataset_df = read_tsv(file='data/TPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv', 
                      col_types = 'fcfffiddiciccccc')

top_cutoff = c(20)
disease_network = lapply(top_cutoff, function(top){
  tmp_df = read_tsv(file = paste0('results/DiseasePathwayNetwork/DiseasePathwayNetwork-FDR<5e-2_top',top,'GOaverageBMA.tsv.gz'), col_types = cols())
  tmp_df$TOP = top
  return(tmp_df)
}) %>% data.table::rbindlist()

disease_pathways = unlist(
  apply(dataset_df,1,function(x){
    c( x[['TargetPathway']],
       disease_network %>% 
         filter(A_name == x[['TargetPathway']]) %>% 
         pull(B_name)
    )
  }) 
)
disease_pathways_freq = data.frame(table(disease_pathways)/nrow(dataset_df)) %>% 
  arrange(-Freq) %>% 
  rename(pathway_name = disease_pathways) %>%
  mutate(Type = ifelse(pathway_name %in% dataset_df$TargetPathway, 'Target','Target-related'))

disease_pathways_freq$pathway_name = factor(disease_pathways_freq$pathway_name,levels = disease_pathways_freq$pathway_name)

ggplot(disease_pathways_freq, aes(x=pathway_name, y=Freq, fill=Type)) +
  geom_bar(position="dodge", stat="identity",width = 0.6, color = "black")+
  xlab('')+
  ylab('Benchmarked Pathway Frequency') +
  scale_fill_manual(values=list(Target='#cc2936',`Target-related`='#08415c')) +
  # ggtitle(paste0(length(unique(disease_pathways_freq$pathway_name)),'/318 pathways in TP benchmark')) + 
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(),axis.line=element_line(colour='black'),
        panel.spacing=unit(1, 'lines'),
        plot.title=element_text(size=13, face='bold'),
        axis.title.y=element_text(size=13,face='bold',colour="black"),
        axis.text.y=element_text(size=10,face='bold',colour="black"),
        axis.text.x=element_text(size=8, face='bold', angle=60, hjust =1,colour="black"),
        strip.text = element_text(size = 12, face = 'bold'),
        legend.justification = c("center", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.position = "top", 
        legend.box = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size=10,face='bold',colour="black"),
        plot.margin=unit(c(1,1,1,1), 'cm'))



##############################
## Measure disease class similarity
##############################

## Subset disease of benchmark
dataset_df = read_tsv(file='data/TPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv', 
                      col_types = 'fcfffiddiciccccc') %>% 
  select(GEO,TargetPathway,) %>%
  droplevels()

top_cutoff = c(10,20,40,'ALL')
pathway_network_df = lapply(top_cutoff, function(top){
  tmp_df = read_tsv(file = paste0('results/DiseasePathwayNetwork/DiseasePathwayNetwork-FDR<5e-2_top',top,'GOaverageBMA.tsv.gz'), col_types = cols())
  
  tmp_df$TOP = top
  return(tmp_df)
}) %>% data.table::rbindlist()
pathway_network_df = pathway_network_df %>%
  filter(TOP == '20') %>% droplevels()

## create combinations
tmp_pathway_combinations = 
  lapply(
    utils::combn(unique(pathway_network_df$A_name), 2, simplify = F),
    function(x){
      return(data.frame(P1 = x[[1]],
                        P2 = x[[2]]))
    }) %>% data.table::rbindlist()

compute_similarity = function(i){
  library(tidyverse)
  x = tmp_pathway_combinations[i,]
  p1 = x[[1]]
  p2 = x[[2]]
  p1_l = pathway_network_df %>% 
    filter(A_name==p1) %>% 
    pull(B_name) %>% 
    unique()
  p2_l = pathway_network_df %>% 
    filter(A_name==p2) %>% 
    pull(B_name) %>% 
    unique()
  ji = length(intersect(p1_l,p2_l)) / length(union(p1_l,p2_l))
  return(
    data.frame(P1=p1,P2=p2,ji=ji)
  )
}

pb = utils::txtProgressBar(min=0, max=nrow(tmp_pathway_combinations), style = 3)
progress = function(n) utils::setTxtProgressBar(pb, n)
opts = list(progress = progress)
cl = parallel::makeCluster(10)
doSNOW::registerDoSNOW(cl)
boot = foreach::foreach(i = 1:nrow(tmp_pathway_combinations), .options.snow = opts)
disease_similarity = foreach::`%dopar%`(boot, compute_similarity(i)) %>% 
  data.table::rbindlist()
parallel::stopCluster(cl)

write_tsv(disease_similarity,'results/DiseasePathwayNetwork/DiseasePathwayNetwork-FDR<5e-2_TOP20-DiseaseSimilarity.tsv.gz')


################# Some stats within subclass
tmp_df = read_tsv(paste0('results/DiseasePathwayNetwork/DiseasePathwayNetwork-FDR<5e-2_TOP20-DiseaseSimilarity.tsv.gz')) 
kegg_df = read_tsv(file = 'data/KEGG/keggrest_hsa.tsv', col_types = 'ffcccff') %>%
  select(pathway_name,subclass) %>% 
  unique() %>%
  droplevels()

tmp_df = tmp_df %>% 
  rename(pathway_name=P1) %>% 
  left_join(kegg_df, by = 'pathway_name') %>% 
  rename(P1=pathway_name,
         pathway_name=P2,
         subclassA=subclass) %>% 
  left_join(kegg_df, by = 'pathway_name') %>% 
  rename(P2=pathway_name,
         subclassB=subclass)

tmp_df %>% 
  filter(subclassA==subclassB) %>% 
  group_by(subclassA) %>% 
  summarise(meanJI = round(mean(ji),2)) %>% 
  arrange(-meanJI) %>% 
  write_tsv(file = 'notebooks/tables/DiseasePathwayNetwork_TOP20-similarity.tsv')
