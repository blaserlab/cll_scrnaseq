# libraries-------------------------------------
library(blaseRtools)# devtools::install_github("git@github.com:blaserlab/blaseRtools.git")
library(tidyverse)
library(monocle3)
library(circlize)
library(ComplexHeatmap)
library(lazyData)
library(cowplot)
library(RColorBrewer)
library(fastSave)
library(ggrepel)

# project data-------------------------------------
# run once to load, run again to unload
requireData(package = "cll.scrnaseq.datapkg",quietly = F)
