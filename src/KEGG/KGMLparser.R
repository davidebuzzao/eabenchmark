#!/usr/bin/env Rscript

########################################################################
## KGML parser
########################################################################
rm(list = ls())
setwd('~/eabenchmark/') 

#curl -s http://rest.kegg.jp/list/pathway/hsa | awk '{split($1','a,":"); print "curl -s http://rest.kegg.jp/get/"a[2]"/kgml -o data/tldkegg/xml/"a[2]".xml"}' | bash

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(parallel))

option_list = list(
  make_option(c("-c", "--ncores"), type="integer", default='8',
              help=" (by default: %default)", metavar="character"),
  make_option(c("-d", "--database"), type="character", default='metabolic',
              help=" (by default: %default)", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
ncores = opt$ncores
db = opt$database

species_id = c('dre','bta','cfa','bsu','mmu','hsa','osa','ssc',
               'spo','dme','sce','rno','ddi','eco','mja',
               'sso','cin','mtu','gga','ath','pfa','cel')

extract_kegg = function(s, db = 'metabolic'){
  library(KEGGgraph)
  library(tidyverse)
  
  if (db == 'metabolic'){
    ###########################################################
    ## METABOLIC
    pathway_id = 
               c('01200','01210','01212','01230','01220','00010','00020','00030',
                 '00040','00051','00052','00053','00500','00520','00620','00630',
                 '00640','00650','00660','00562','00190','00195','00196','00710',
                 '00720','00680','00910','00920','00061','00062','00071','00072',
                 '00073','00100','00120','00121','00140','00561','00564','00565',
                 '00600','00590','00591','00592','00230','00240','00250','00260',
                 '00270','00280','00290','00300','00310','00220','00330','00340',
                 '00350','00360','00380','00400','00410','00430','00440','00450',
                 '00460','00471','00472','00473','00480','00510','00513','00512',
                 '00515','00514','00532','00534','00533','00531','00563','00601',
                 '00603','00604','00540','00550','00511','00571','00572','00730',
                 '00740','00750','00760','00770','00780','00785','00790','00670',
                 '00830','00860','00130','00900','00902','00909','00904','00906',
                 '00905','00981','00908','00903','00281','00522','00253','00523',
                 '00940','00945','00941','00944','00942','00943','00901','00403',
                 '00950','00960','00232','00965','00966','00402','00311','00332',
                 '00261','00331','00521','00524','00525','00401','00404','00405',
                 '00333','00254','00999','00998','00997','00362','00627','00364',
                 '00625','00361','00623','00622','00633','00642','00643','00791',
                 '00930','00363','00621','00626','00624','00365','00984','00980',
                 '00982','00983','00970','01040','01051','01052','01053','01054',
                 '01055','01056','01057','01058','02010','02060','03008','03010',
                 '03013','03015','03018','03020','03022','03030','03040','03050',
                 '03060','03070','03410','03420','03430','03440','03450','03460',
                 '04120','04122','04130','04139','04140','04141','04142','04144',
                 '04145','04146','04973','04974','04975','04977','07222','01232',
                 '01250','01240','00470','00542','00541','00552','00543','00946','00996')
    
  } else if (db == 'signaling'){
    ###########################################################
    ## SIGNALING
    pathway_id =
               c('04130','04014','04015','04010','04013','04016','04011','04012',
                 '04310','04330','04340','04341','04350','04390','04391','04392',
                 '04370','04371','04630','04064','04668','04066','04068','04020',
                 '04070','04072','04071','04024','04022','04151','04152','04150',
                 '04075','04080','04060','04061','04512','04514','04115','04510',
                 '04520','04530','04540','02024','04620','04624','04621','04622',
                 '04623','04625','04660','04657','04662','04664','04062','04910',
                 '04922','04920','03320','04912','04915','04917','04921','04926',
                 '04919','04260','04261','04270','04720','04723','04722','04744',
                 '04745','04740','04742','04626','02020','04110','04111','04112',
                 '04113','04114','04210','04610','04612','04614','04640','04650',
                 '04666','04670','04672','04721','04724','04725','04726','04727',
                 '04728','04730','04914','04916')
  }
  
  message(s)
  dirout = paste0('data/kegg/xml/',s)
  if(!dir.exists(dirout)){
    dir.create(dirout)
    print(paste("Folder created at", dirout))
  } else {
    print(paste("Folder already exists at", dirout))
  }
  
  f1out = paste0(dirout,'/', s,'.tsv')
  if(!file.exists(f1out)){
    message("Fetching pathway list")
    system(paste0('wget --quiet --output-document ',f1out, ' http://rest.kegg.jp/list/pathway/',s))
  }
  tmp_df = read_tsv(file = f1out, col_names = c('pathway_id','pathway_name'), col_types = 'cc')
  tmp_df$pathway_id = gsub(paste0('path:',s),'', tmp_df$pathway_id)
  tmp_pathway_id = pathway_id[pathway_id %in% tmp_df$pathway_id]
  
  count = 1
  tot = length(tmp_pathway_id)
  print("Downloading pathways")
  for (p in tmp_pathway_id){
    message(paste(count,tot,sep = '/')); count = count + 1
    f2out = paste0(dirout,'/',s,p,'.xml')
    if(!file.exists(f2out)){
      system(paste0('wget --quiet --output-document ',f2out, ' http://rest.kegg.jp/get/',s,p,'/kgml'))
    }
  }

  type_rel = c('PPrel','ECrel','GErel')
  
  count = 1
  tot = length(tmp_pathway_id)
  message("Building pathway networks")
  df = lapply(tmp_pathway_id,function(x){
    print(paste(count,tot,sep = '/')); count <<- count + 1
    # message(x[['pathway_id']])
    pathway_file = paste0('data/kegg/xml/',s,'/',s,x,'.xml')
    if (!file.exists(pathway_file)) return()
    pathway_network = 
      tryCatch(
        parseKGML2DataFrame(pathway_file, reactions = F),
        error=function(con){
          return(NULL)})
    if(nrow(pathway_network)==0 | is.null(pathway_network)) return()
    
    pathway_network$pathway_id = x
    pathway_network$from = gsub(paste0(s,':'),'',pathway_network$from)
    pathway_network$to = gsub(paste0(s,':'),'',pathway_network$to)
    pathway_network = pathway_network %>% filter(type %in% type_rel)
    if(nrow(pathway_network)==0 | is.null(pathway_network)) return()
    
    pathway_network = apply(pathway_network,1,function(x){
      # print(x)
      p = sort(c(x[['from']],x[['to']]))
      x[['from']] = p[[1]]
      x[['to']] = p[[2]]
      return(data.frame(p1 = p[[1]],
                        p2 = p[[2]],
                        type = x[['type']],
                        subtype = x[['subtype']],
                        pathway_id = x[['pathway_id']]))
    }) %>% data.table::rbindlist()
    
    return(pathway_network)
  }) %>% data.table::rbindlist() 
  
  if (nrow(df)==0){
    return(
      data.frame(Species = s,
                 Genes = 0,
                 Links = 0,
                 Pathways = 0)
    )
  } else {
    k = length(unique(df$pathway_id))
    g = length(unique(c(df$p1,df$p2)))
    
    df = df %>% 
      filter(p1!=p2) %>%
      group_by(p1,p2) %>% 
      summarise(pathway_id = paste(pathway_id, collapse = ','),
      )
    l = nrow(df)
    
    df %>% 
      filter(p1!=p2) %>% 
      write_tsv(paste0('data/kegg/xml/',s,'_KEGGgraph.tsv'))
    
    return(
      data.frame(Species = s,
                 Genes = g,
                 Links = l,
                 Pathways = k)
    )
  }
}


#################################################
## EXECUTE METHOD
pb = utils::txtProgressBar(min=0, max=length(species_id), style = 3)
progress = function(n) utils::setTxtProgressBar(pb, n)
opts = list(progress = progress)
cl = parallel::makeCluster(ncores)
doSNOW::registerDoSNOW(cl)
boot = foreach::foreach(i = species_id, .options.snow = opts)
kegg_df = foreach::`%dopar%`(boot, extract_kegg(i,db)) %>% 
  data.table::rbindlist()
parallel::stopCluster(cl)

kegg_df %>% 
  write_tsv(paste0('data/KEGG/xml/KEGG-',db,'.tsv'))
