# source("R/dependencies.R")
# source("R/configs.R")
source("R/fig_staging/umap_partition_assignment.R")
source("R/fig_staging/umap_density.R")
source("R/fig_staging/cell_representation_barchart.R")
source("R/fig_staging/umap_da_score.R")
source("R/fig_staging/umap_leiden_enrichment.R")
source("R/fig_staging/leiden_enrichment_tm_hm.R")

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
    umap_da_score,
    umap_leiden_enrichment,
    labels = c("C","D", "E"),
    align = "h",
    axis = "b",
    ncol = 3,
    rel_widths = c(0.67, 1, 1)
    )

leiden_enrichment_tm_hm_alt <- plot_grid(leiden_enrichment_tm_hm) + theme(plot.margin = unit(c(0,5,0,0), "mm"))    
fig_1_mid2 <- 
  plot_grid(
    leiden_enrichment_tm_hm_alt,
    labels = c("F"),
    align = "h",
    axis = "b",
    ncol = 1
  )


fig_1_bottom <- 
  plot_grid(
    NULL,
    NULL,
    ncol = 2,
    rel_widths = c(1,1),
    labels = c("G","H")
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
 filename = "test.png",
 # filename = fs::path(network_out, "fig_1.png"),
 base_width = 7.5,
 base_height = 9.75
)
