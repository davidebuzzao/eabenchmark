
run_binox = function(GEO,
                     network_name,
                     benchmark='TPbenchmark',
                     binox_dir){
  indir = getwd()
  network_file = paste0(indir, 'data/binox/BinoX-',network_name,'.tsv.randNet')
  
  #################################################
  ## EXECUTE METHOD
  system(paste0('/usr/bin/time -f \'%E %U %S\' -o results/binox_runtime.txt bash -c \'',
                binox_dir,'/BinoX -r ',network_file,' ',
                '-a <(tail -n +2 ', indir, 'data/methodSpecific/binox/kegg-',network_name,'.tsv | cut -f 1,2) ',
                '-b <(tail -n +2 ', indir, 'data/methodSpecific/binox/',benchmark,'/',network_name,'/',GEO,'.tsv | cut -f 1,2) ',
                '-o ', indir, 'results/',GEO,'_binox-',network_name,'.tsv ',
                '-p large ',
                '-g 1\''))
  
  enrichment_df = read_tsv(file = paste0('results/',GEO,'_binox-',network_name,'.tsv'))
  return(enrichment_df)
}

