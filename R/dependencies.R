# set up the renv and repair with snapshot if needed

# blaseRtemplates::easy_install("ggrastr", "link_from_cache")

# load core packages for the analysis
suppressPackageStartupMessages(library("conflicted"))
suppressPackageStartupMessages(library("blaseRtools"))
suppressPackageStartupMessages(library("blaseRdata"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("monocle3"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("lazyData"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("rstatix"))
suppressPackageStartupMessages(library("pander"))
suppressPackageStartupMessages(library("Rcpp"))
suppressPackageStartupMessages(library("ggrastr"))

# uncomment and use the following to install or update the data package---------------------------------------
blaseRtemplates::bb_renv_datapkg(path = "~/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/datapkg")

# project data-------------------------------------
# run once to load, run again to unload
requireData(package = "cll.scrnaseq.datapkg",quietly = F)
