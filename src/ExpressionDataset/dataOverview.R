rm(list = ls())
setwd('~/eabenchmark/') 

library(tidyverse)
dataset_df = read_tsv(file='data/TPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv', 
                      col_types = 'fcfffiddiciccccc')

plot_overview =  function(row_n){
  library(tidyverse)
  library(SummarizedExperiment)
  library(gridExtra)
  library(grid)
  library(ggrepel)
  
  x = dataset_df[row_n,]
  alpha = as.double(x[['alpha']])
  beta = as.double(x[['beta']])
  num_deg = x[['NrOfDEG']]
  
  se = readRDS(paste0('data/TPbenchmark/',x[['DataType']],'/',gsub('\\.[0-9]*','',x[['GEO']]),'.Rds'))
  se@colData$GROUP = as.numeric(as.character(se@colData$GROUP))
  colnames(SummarizedExperiment::rowData(se)) = c('FC','p.value','t','ADJ.PVAL')
  se@metadata$annotation = 'hsa'
  se@metadata$dataType = x[['DataType']]

  if(any(is.na(SummarizedExperiment::assay(se)))){
    ## GSE54236 --> it's the only dataset to have NaN values and these methods can't coope with that 
    se = subset(se,complete.cases(assay(se)))
  }
  
  ano = colData(se)
  expression = assay(se)
  
  control_samples = rownames(ano[ano$GROUP==0,])
  disease_samples =  rownames(ano[ano$GROUP==1,])
  
  cc_dataset = cor(expression, method = 'spearman')
  cc_control = cor(expression[,control_samples], method = 'spearman')
  cc_disease = cor(expression[,disease_samples], method = 'spearman')
  cc_dataset.df = stack(cc_dataset)
  cc_control.df = stack(cc_control)
  cc_disease.df = stack(cc_disease)
  
  cc_dataset.df$Dataset = 'All'
  cc_control.df$Dataset = 'Control'
  cc_disease.df$Dataset = 'Disease'
  cc.df = as.data.frame(rbind(cc_dataset.df,cc_control.df,cc_disease.df))
  
  dp = ggplot(cc.df, aes(x = value, color=Dataset)) +
    geom_density() + 
    xlab('Spearman correlation (Expression)') + 
    scale_color_manual(values = c(
      'All' = '#277da1',
      'Control' = '#a7c957',
      'Disease' = '#c1121f'
    )) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = 'black'),
          plot.title = element_text(size=12, face='bold'),
          axis.text.y = element_text(size=12, face='bold'),
          axis.text.x  = element_text(size=12, face='bold'),
          axis.title = element_text(size=12, face='bold'),
          legend.title = element_text(size=12, face='bold'),
          legend.text = element_text(size=12),
          legend.justification = c('right', 'top'),
          legend.box.just = 'right',
          legend.margin = margin(6, 6, 6, 6))

  dat_bp = stack(as.data.frame(expression))
  dat_bp$Group = as.factor(ifelse(dat_bp$ind %in% control_samples, 'Control', 'Disease'))
  bp = ggplot(dat_bp, aes(fill=Group)) +
    geom_boxplot(aes(x = ind, y = values), outlier.size = 0.5) +
    scale_fill_manual(values = c(
      'Control' = '#a7c957',
      'Disease' = '#c1121f'
    )) +
    xlab('Samples') + ylab('Expression') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = 'black'),
          plot.title = element_text(size=12, face='bold'),
          axis.text.y = element_text(size=12, face='bold'),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size=12, face='bold'),
          legend.title = element_text(size=12, face='bold'),
          legend.text = element_text(size=12),
          legend.justification = c('right', 'top'),
          legend.box.just = 'right',
          legend.margin = margin(6, 6, 6, 6))

  #Sample PCA
  dat_pca = t(expression)
  
  pcaRes = tryCatch(
    {prcomp(dat_pca)},
    error=function(cond) {
      return(NULL)})
  
  if (is.null(pcaRes)){
    print('Error in computing PCA')
    next}
  varExp = round(pcaRes$sdev^2 / sum(pcaRes$sdev^2) * 100)
  pcaDF = data.frame(
    PC1 = pcaRes$x[, 1],
    PC2 = pcaRes$x[, 2],
    PC3 = pcaRes$x[, 3],
    Group = as.factor(ifelse(ano$GROUP==0,'Control','Disease')),
    Sample = rownames(ano))
  
  pc12 = ggplot(
    data = pcaDF,
    mapping = aes(x = PC1, y = PC2, color = Group, label = Sample)) +
    geom_point(size = 3) +
    scale_color_manual(values = c(
      'Control' = '#a7c957',
      'Disease' = '#c1121f'
    )) +
    labs(x = paste0('PC1 (', varExp[1], ' %)'), y = paste0('PC2 (', varExp[2], ' %)')) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = 'black'),
          plot.title = element_text(size=12, face='bold'),
          axis.text.y = element_text(size=12, face='bold'),
          axis.text.x  = element_text(size=12, face='bold'),
          axis.title = element_text(size=12, face='bold'),
          legend.title = element_text(size=12, face='bold'),
          legend.text = element_text(size=12),
          legend.justification = c('right', 'top'),
          legend.box.just = 'right',
          legend.margin = margin(6, 6, 6, 6))
  
  pc23 = ggplot(
    data = pcaDF, mapping = aes(x = PC2, y = PC3, color = Group, label = Sample)) +
    geom_point(size = 3) +
    scale_color_manual(values = c(
      'Control' = '#a7c957',
      'Disease' = '#c1121f'
    )) +
    labs(x = paste0('PC2 (', varExp[2], ' %)'), y = paste0('PC3 (', varExp[3], ' %)')) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = 'black'),
          plot.title = element_text(size=12, face='bold'),
          axis.text.y = element_text(size=12, face='bold'),
          axis.text.x  = element_text(size=12, face='bold'),
          axis.title = element_text(size=12, face='bold'),
          legend.title = element_text(size=12, face='bold'),
          legend.text = element_text(size=12),
          legend.justification = c('right', 'top'),
          legend.box.just = 'right',
          legend.margin = margin(6, 6, 6, 6))
  
  ##################################
  ### Differential expression analysis
  dea = as.data.frame(rowData(se))
  # Classify probes for plotting
  dea$diffexpr = 'No'
  dea$diffexpr[dea$ADJ.PVAL < alpha & sign(dea$t)==1] = 'Up'
  dea$diffexpr[dea$ADJ.PVAL < alpha & sign(dea$t)==-1] = 'Down'
  # Add gene symbols as labels to DE genes
  dea$delabel = rownames(dea)
  dea$delabel[dea$diffexpr == 'No'] = ''

  # Volcano plot for bipolar patients vs controls
  vp1 = ggplot(data = as.data.frame(dea),aes(x = FC,y = -log10(ADJ.PVAL),color = diffexpr,label = delabel)) +
    geom_point() +
    geom_hline(yintercept = -log10(alpha), col = 'gray45', linetype = 'dashed') +
    labs(x = 'log2(FoldChange)', y = '-log10(FDR)') +
    scale_color_manual(values = c('Down'='#277da1', 
                                  'No'='#353535', 
                                  'Up'='#c1121f')) +
    geom_text_repel(show.legend = FALSE) +
    ggtitle('~Group') + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = 'black'),
          plot.title = element_text(size=12, face='bold'),
          axis.text.y = element_text(size=12, face='bold'),
          axis.text.x  = element_text(size=12, face='bold'),
          axis.title = element_text(size=12, face='bold'),
          legend.title = element_text(size=12, face='bold'),
          legend.text = element_text(size=12),
          legend.justification = c('right', 'top'),
          legend.box.just = 'right',
          legend.margin = margin(6, 6, 6, 6))      
    
    g = grid.arrange(bp,dp,
                     pc12,pc23,vp1,
                     nrow = 2,
                     top = textGrob(paste0(x[['GEO']],' (',x[['TargetPathway']],')'),gp=gpar(fontsize=20,font=3)),
                     layout_matrix = cbind(c(1,3), c(1,4), c(2,5)))
    
    print(g)
    ggsave(paste0('data/TPbenchmark/statsOverview/',x[['GEO']],'.png'),  plot = g,
           width = 16, height = 8, dpi = 300, units = "in", device='png')
 }

#################################################
## EXECUTE METHOD
pb = utils::txtProgressBar(min=0, max=nrow(dataset_df), style = 3)
progress = function(n) utils::setTxtProgressBar(pb, n)
opts = list(progress = progress)
cl = parallel::makeCluster(8)
doSNOW::registerDoSNOW(cl)
boot = foreach::foreach(i = 1:nrow(dataset_df), .options.snow = opts)
foreach::`%dopar%`(boot, plot_overview(i)) 
parallel::stopCluster(cl)

