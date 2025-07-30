# source("R/dependencies.R")
# source("R/configs.R")

source_all_figs <- TRUE
# source_all_figs <- FALSE

if (source_all_figs) {
  source("R/fig_staging/mono_seurat_celltypes.R")
  source("R/fig_staging/mono_louvain_umap.R")
  source("R/fig_staging/mono_cluster_assignment.R")
  source("R/fig_staging/cd14_enrichment_umap.R") 
  source("R/fig_staging/cd16_enrichment_umap.R") 
  
}

fig_S7_top <- plot_grid(
  mono_seurat_celltypes + theme(aspect.ratio = 0.85),
  mono_louvain_umap + theme(aspect.ratio = 0.85),
  ncol = 2, 
  align = "h", 
  axis = "b",
  labels = c("A", "B"),
  rel_widths = c(1, 1)
)
fig_S7_mid <- plot_grid(
  mono_cluster_assignment,
  NULL,
  align = "h",
  axis = "b",
  labels = c("C", ""),
  ncol = 2,
  rel_widths = c(1,0.5)
)

fig_S7_bot <- plot_grid(
  cd14_enrichment_umap + theme(aspect.ratio = 0.9),
  cd16_enrichment_umap + theme(aspect.ratio = 0.9),
  labels = c("D", "E"),
  ncol =  2,
  rel_widths = c(1,1)
)

fig_S7 <- plot_grid(
  fig_S7_top,
  fig_S7_mid,
  fig_S7_bot,
  ncol = 1,
  rel_heights = c(1,1,1)
)

save_plot(
  fig_S7,
  filename = fs::path(network_out, "fig_S7.pdf"),
  base_width = 7.5,
  base_height = 9.75
)
