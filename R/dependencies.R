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
suppressPackageStartupMessages(library("tidymodels"))
suppressPackageStartupMessages(library("nlme"))
suppressPackageStartupMessages(library("multilevelmod"))
suppressPackageStartupMessages(library("patchwork"))



# uncomment and use the following to install or update the data package---------------------------------------
blaseRtemplates::project_data(path = "~/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/datapkg")

