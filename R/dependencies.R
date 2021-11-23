# set up the renv and repair with snapshot if needed
# renv::init()
# renv::snapshot()

# uncomment and run to restore packages from the lock file (for example in new installation)
#renv::restore()

# blaseRtools and additional dependencies you may have to install since they are not recognized by renv::init
# renv::install("/usr/lib/R/site-library/blaseRtools")

# load core packages for the analysis
library("blaseRtools")
library("tidyverse")
library("monocle3")
library("circlize")
library("ComplexHeatmap")
library("lazyData")
library("cowplot")
library("RColorBrewer")
library("ggrepel")
library("ggpubr")
library("rstatix")
library("pander")
library("Rcpp")
# uncomment and use the following to install or update the data package---------------------------------------
bb_renv_datapkg(path = "~/network/X/Labs/Blaser/cll_scrnaseq_manuscript/datapkg")

# project data-------------------------------------
# run once to load, run again to unload
requireData(package = "cll.scrnaseq.datapkg",quietly = F)
