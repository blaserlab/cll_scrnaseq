# set up the renv and repair with snapshot if needed
# renv::init()
# renv::snapshot()

# uncomment and run to restore packages from the lock file (for example in new installation)
#renv::restore()

# blaseRtools and additional dependencies you may have to install since they are not recognized by renv::init
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/blaseRtools")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/DESeq2")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/genefilter")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/annotate")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/AnnotationDbi")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/KEGGREST")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/Biostrings")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/geneplotter")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/DoubletFinder")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/Seurat")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/SeuratDisk")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/rrvgo")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/GO.db")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/GOSemSim")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/scater")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/topGO")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/fastSave")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/hexbin")
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/graph")

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

# uncomment and use the following to install or update the data package---------------------------------------
# bb_renv_datapkg(path = "~/network/X/Labs/Blaser/cll_scrnaseq_manuscript/datapkg")

# project data-------------------------------------
# run once to load, run again to unload
requireData(package = "cll.scrnaseq.datapkg",quietly = F)
