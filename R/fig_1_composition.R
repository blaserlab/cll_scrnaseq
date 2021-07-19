source("R/configs.R")
source("R/fig_1_staging.R")

fig_1_top <-
  plot_grid(
    umap_density,
    NULL,
    ncol = 2,
    rel_widths = c(4,1),
    labels = c("A","")
    )
    
fig_1_mid <- 
  plot_grid(
    volcano_BTK_bcells,
    subpop_top_markers_heatmap,
    labels = c("B","C"),
    align = "h",
    axis = "b",
    ncol = 2,
    rel_widths = c(1,2)
  )



fig_1_bottom <- 
  plot_grid(
    cluster_proportion_ratio_plot,
    module_heatmap_bcells,
    ncol = 2,
    rel_widths = c(1,1),
    labels = c("D","E")
  )

fig_1 <- 
  plot_grid(
    fig_1_top,
    fig_1_mid,
    fig_1_bottom,
    nrow = 3,
    rel_heights = c(1.2,1,1.5)
  )

save_plot(
 plot = fig_1,
 # filename = "test.png",
 filename = str_glue("{network_out}/fig_1.png"),
 base_width = 7.5,
 base_height = 9.75
)