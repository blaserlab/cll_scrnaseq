# source("R/configs.R")
# source("R/fig_1_staging.R")

fig_1_top <-
  plot_grid(
    umap_partition_assignment,
    NULL,
    ncol = 2,
    rel_widths = c(4,1),
    labels = c("A","")
    )

fig_1_mid1 <- 
  plot_grid(
    volcano_BTK_bcells,
    umap_subcluster,
    NULL,
    labels = c("B","C", ""),
    align = "h",
    axis = "b",
    ncol = 3
  )

subpop_top_markers_heatmap_alt <- plot_grid(subpop_top_markers_heatmap) + theme(plot.margin = unit(c(0,5,0,0), "mm"))    
fig_1_mid2 <- 
  plot_grid(
    subpop_top_markers_heatmap_alt,
    labels = c("D"),
    align = "h",
    axis = "b",
    ncol = 1
  )


module_heatmap_bcells_alt <- plot_grid(module_heatmap_bcells) + theme(plot.margin = unit(c(0, 5, 0, 0), "mm"))
fig_1_bottom <- 
  plot_grid(
    cluster_proportion_ratio_plot,
    module_heatmap_bcells_alt,
    ncol = 2,
    rel_widths = c(1,1),
    labels = c("E","F")
  )

fig_1 <- 
  plot_grid(
    fig_1_top,
    fig_1_mid1,
    fig_1_mid2,
    fig_1_bottom,
    nrow = 4,
    rel_heights = c(1.5,1,0.8,1.2)
  )

save_plot(
 plot = fig_1,
 # filename = "test.png",
 filename = fs::path(network_out, "fig_1.png"),
 base_width = 7.5,
 base_height = 9.75
)