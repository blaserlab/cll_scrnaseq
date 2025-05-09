# source("R/dependencies.R")
# source("R/configs.R")

# source_all_figs <- TRUE
source_all_figs <- FALSE

if (source_all_figs) {
  source("R/fig_staging/tnkcell_subpop_umap.R")
  source("R/fig_staging/umap_da_score_tnk.R")
  source("R/fig_staging/cluster_proportion_ratio_plot_tnk2.R")
  source("R/fig_staging/abs_treg.R")
  source("R/fig_staging/tcell_gsea.R")
  source("R/fig_staging/tcell_leiden_enrichment_genebub.R")
  
}

fig_3_top <- plot_grid(
  tnkcell_subpop_umap + theme(aspect.ratio = 0.85),
  umap_da_score_tnk + theme(aspect.ratio = 0.85),
  ncol = 2, 
  align = "h", 
  axis = "b",
  labels = c("A", "B"),
  rel_widths = c(1.25, 1)
)
fig_3_mid <- plot_grid(
  cluster_proportion_plot_tnk2,
  tcell_gsea_plots_combined,
  
  NULL,
  align = "h",
  axis = "b",
  labels = c("C", "D", ""),
  ncol = 3,
  rel_widths = c(1,1,0.25)
)

fig_3_bot <- plot_grid(
  tcell_leiden_enrichment_genebub,
  abs_treg,
  align = "h",
  axis = "b",
  labels = c("E", "F"),
  ncol =  2,
  rel_widths = c(1, 1)
)

fig_3 <- plot_grid(
  fig_3_top,
  fig_3_mid,
  fig_3_bot,
  ncol = 1,
  rel_heights = c(1,1,1)
)

save_plot(
  fig_3,
  filename = fs::path(network_out, "fig_3.pdf"),
  base_width = 7.5,
  base_height = 9.75
)
