# uncomment and install if necessary
# install.packages('blaseRdata', repos = c('https://blaserlab.r-universe.dev', 'https://cloud.r-project.org'))
# install.packages('blaseRtools', repos = c('https://blaserlab.r-universe.dev', 'https://cloud.r-project.org'))

# load core packages for the analysis
suppressPackageStartupMessages(library("conflicted"))
suppressPackageStartupMessages(library("blaseRtools"))
suppressPackageStartupMessages(library("blaseRdata"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("monocle3"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("rstatix"))
suppressPackageStartupMessages(library("pander"))
suppressPackageStartupMessages(library("Rcpp"))
suppressPackageStartupMessages(library("ggrastr"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library(lmerTest))  # lmer + Satterthwaite/Kenward-Roger p-values
suppressPackageStartupMessages(library(emmeans))  # estimated marginal means & contrasts


project_data(path = "~/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/datapkg")

