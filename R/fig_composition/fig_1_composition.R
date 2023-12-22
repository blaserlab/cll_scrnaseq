# source("R/dependencies.R")
# source("R/configs.R")

source_all_figs <- TRUE 
# source_all_figs <- FALSE

if (source_all_figs) {
  source("R/fig_staging/umap_partition_assignment.R")
  source("R/fig_staging/umap_density.R")
  source("R/fig_staging/cell_representation_barchart.R")
  source("R/fig_staging/umap_da_score.R")
  source("R/fig_staging/umap_leiden_enrichment.R")
  source("R/fig_staging/leiden_enrichment_tm_hm.R")
  source("R/fig_staging/umap_numbered_leiden.R")
  source("R/fig_staging/cluster_proportion_ratio_plot.R")
  source("R/fig_staging/umap_subcluster.R")
  
}

fig_1_top <-
  plot_grid(
    umap_partition_assignment,
    umap_density,
    ncol = 2,
    rel_widths = c(0.67,1),
    labels = c("A","B")
    )

fig_1_mid1 <- 
  plot_grid(
    cell_representation_barchart,
    umap_numbered_leiden,
    umap_da_score,
    labels = c("C","D", "E"),
    align = "h",
    axis = "b",
    ncol = 3,
    rel_widths = c(0.8, 1, 0.8)
    )


fig_1_mid2 <- 
  plot_grid(
    umap_leiden_enrichment,
    umap_subcluster,
    cluster_proportion_ratio_plot,
    labels = c("F", "G", "H"),
    align = "h",
    axis = "b",
    ncol = 3,
    rel_widths = c(0.8, 1, 0.8)
  )


leiden_enrichment_tm_hm_alt <- plot_grid(leiden_enrichment_tm_hm) + theme(plot.margin = unit(c(0,5,0,0), "mm"))    
fig_1_bottom <- 
  plot_grid(
    leiden_enrichment_tm_hm_alt,
    ncol = 1,
    labels = c("I")
  )

fig_1 <- 
  plot_grid(
    fig_1_top,
    fig_1_mid1,
    fig_1_mid2,
    fig_1_bottom,
    nrow = 4,
    rel_heights = c(1,1,1,1)
  )

save_plot(
 plot = fig_1,
 # filename = "fig_1.png",
 filename = fs::path(network_out, "fig_1.pdf"),
 base_width = 7.5,
 base_height = 9.75
)
