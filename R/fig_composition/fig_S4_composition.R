# source("R/dependencies.R")
# source("R/configs.R")

source_all_figs <- TRUE
# source_all_figs <- FALSE

if (source_all_figs) {
  source("R/fig_staging/tnk_seurat_celltypes.R")
  source("R/fig_staging/tnk_leiden_umap.R")
  source("R/fig_staging/tnk_cluster_assignment.R")
  source("R/fig_staging/umap_subcluster_tnk.R")
  
}

fig_S4_top <- plot_grid(
  tnk_seurat_celltypes + theme(aspect.ratio = 0.85),
  tnk_leiden_umap + theme(aspect.ratio = 0.85),
  ncol = 2, 
  align = "h", 
  axis = "b",
  labels = c("A", "B"),
  rel_widths = c(1, 1)
)
fig_S4_mid <- plot_grid(
  tnk_cluster_assignment,
  NULL,
  align = "h",
  axis = "b",
  labels = c("C", ""),
  ncol = 2,
  rel_widths = c(1,0.5)
)

fig_S4_bot <- plot_grid(
  umap_subcluster_tnk + theme(aspect.ratio = 0.85),
  NULL,
  labels = c("D", ""),
  ncol =  2,
  rel_widths = c(1,1)
)

fig_S4 <- plot_grid(
  fig_S4_top,
  fig_S4_mid,
  fig_S4_bot,
  ncol = 1,
  rel_heights = c(1,1,1)
)

save_plot(
  fig_S4,
  filename = fs::path(network_out, "fig_S4.pdf"),
  base_width = 7.5,
  base_height = 9.75
)