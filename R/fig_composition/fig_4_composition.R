# source("R/dependencies.R")
# source("R/configs.R")

source_all_figs <- TRUE
# source_all_figs <- FALSE

if (source_all_figs) {
  source("R/fig_staging/mono_subpop_umap.R")
  source("R/fig_staging/umap_da_score_mono.R")
  source("R/fig_staging/cluster_proportion_ratio_plot_mono.R")
  source("R/fig_staging/cd14_gsea.R")
  source("R/fig_staging/cd16_gsea.R")
  source("R/fig_staging/mono_genebub.R")
  
}

fig_4_top <- plot_grid(
  mono_subpop_umap + theme(aspect.ratio = 0.85),
  umap_da_score_mono + theme(aspect.ratio = 0.85),
  ncol = 2, 
  align = "h", 
  axis = "b",
  labels = c("A", "B"),
  rel_widths = c(1.25, 1)
)
fig_4_mid <- plot_grid(
  cluster_proportion_plot_mono,
  cd14_gsea_plots_combined,
  NULL,
  align = "h",
  axis = "b",
  labels = c("C", "D", ""),
  ncol = 3,
  rel_widths = c(1,1,0.25)
)

fig_4_bot <- plot_grid(
  cd16_gsea_plots_combined,
  mono_genebub,
  align = "h",
  axis = "b",
  labels = c("E", "F"),
  ncol =  2,
  rel_widths = c(1, 1.5)
)

fig_4 <- plot_grid(
  fig_4_top,
  fig_4_mid,
  fig_4_bot,
  ncol = 1,
  rel_heights = c(1,1,1)
)

save_plot(
  fig_4,
  filename = fs::path(network_out, "fig_4.pdf"),
  base_width = 7.5,
  base_height = 9.75
)
