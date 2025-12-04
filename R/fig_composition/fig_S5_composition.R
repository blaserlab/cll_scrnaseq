make_fig_S5 <- function(source_all_figs = TRUE) {
  if (source_all_figs) {
    source("R/fig_staging/tnk_seurat_celltypes.R", local = TRUE)
    source("R/fig_staging/tnk_leiden_umap.R", local = TRUE)
    source("R/fig_staging/tnk_cluster_assignment.R", local = TRUE)
    source("R/fig_staging/umap_subcluster_tnk.R", local = TRUE)
    source("R/fig_staging/tcr_diversity_plot.R", local = TRUE)
    source("R/fig_staging/cd4_gsea_plots.R", local = TRUE)
  }
  
  fig_S5_top <- plot_grid(
    tnk_seurat_celltypes + theme(aspect.ratio = 0.85),
    tnk_leiden_umap + theme(aspect.ratio = 0.85),
    ncol = 2,
    align = "h",
    axis = "b",
    labels = c("A", "B"),
    rel_widths = c(1, 1)
  )
  fig_S5_mid <- plot_grid(
    tnk_cluster_assignment,
    umap_subcluster_tnk + theme(aspect.ratio = 0.85),
    align = "h",
    axis = "b",
    labels = c("C", "D"),
    ncol = 2,
    rel_widths = c(1, 1)
  )
  
  fig_S5_bot <- plot_grid(
    cd4_tcell_gsea_plots_combined,
    tcr_diversity_plot,
    labels = c("E", "F"),
    ncol =  2,
    rel_widths = c(1, 1)
  )
  
  fig_S5 <- plot_grid(
    fig_S5_top,
    fig_S5_mid,
    fig_S5_bot,
    ncol = 1,
    rel_heights = c(1, 1, 1)
  )
  
  save_plot(
    fig_S5,
    filename = fs::path(network_out, "fig_S5.pdf"),
    base_width = 7.5,
    base_height = 9.75
  )
}
make_fig_S5()