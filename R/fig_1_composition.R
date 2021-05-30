# source("R/configs.R")
# source("fig_1_staging.R")

fig_1_top <-
  plot_grid(
    umap_partitions + theme(plot.margin = margin(t = 15)),
    umap_density,
    align = "h",
    axis = "b",
    ncol = 2, 
    rel_widths = c(1,1),
    labels = c("A","B")
    )
    
fig_1_mid <- 
  plot_grid(
    umap_binned_leiden,
    volcano_MRD1_BTK,
    volcano_MRD2_BTK,
    labels = c("C","D","E"),
    align = "h",
    axis = "b",
    ncol = 3,
    rel_widths = c(1,1,1)
  )



fig_1_bottom <- 
  plot_grid(
    module_heatmap_bcells,
    mod4_bubble + theme(plot.margin = margin(t = 35, b = 35)),
    ncol = 2,
    rel_widths = c(1,1),
    labels = c("E","F")
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
 filename = str_glue("{network_out}/fig_1.png"),
 base_width = 7.5,
 base_height = 9.75
)