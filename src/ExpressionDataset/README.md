# Pipeline to generate benchmark datasets

In this folder you will find the following scripts:
  
  - `select_gemma.R` and `fetch_gemma.R`, to select and fetch data from [Gemma](https://gemma.msl.ubc.ca) database.
  - `dataOverview.R`, to produce summary overview plots of all expression datasets.
  - `topDEG.R`, to select DEG with min 15, max 500 and FDR<0.2.
  - `geneLabels_resampling.R`, to resample genelabels for expression dataset from Uniprot mapping file.
  - `extract_deg.R`, to extract DEG for both TP or FP benchmark.

Note that steps of manual intervention were needed to select the final 82 datasets.