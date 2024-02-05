#!/usr/bin/env Rscript

## CePA 0.8.0
install.packages('CePa',version = "0.8.0",repos = "http://cran.us.r-project.org")

## GSVA v1.42.0
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GSVA")

## NEAT v1.2.3
install.packages('neat', version = '1.2.3',repos = "http://cran.us.r-project.org")

## ANUBIX v1.0.3
require(devtools)
install_bitbucket("sonnhammergroup/anubix@d6ed5b9")
