source("R/configs.R")
source("R/fig_S1_staging.R")

fig_S1_top <- 
  plot_grid(
    umap_partitions,
    umap_bcell_leiden,
    ncol = 2,
    rel_widths = c(1,1),
    align = "h",
    axis = "b",
    labels = c("A","B")
  )

fig_S1_mid <- 
  plot_grid(
    umap_binned_leiden,
    NULL,
    ncol = 2,
    rel_widths = c(1,1),
    labels = c("C","")
  )

fig_S1_bot <- 
  plot_grid(mod4_bubble,
            NULL,
            ncol = 2,
            rel_widths = c(2,1),
            labels = c("D","")
            )


fig_S1 <- 
  plot_grid(
    fig_S1_top,
    fig_S1_mid,
    fig_S1_bot,
    nrow = 3,
    rel_heights = c(1,1,1.5))

save_plot(
 plot = fig_S1,
 filename = str_glue("{network_out}/fig_S1.png"),
 base_width = 7.5,
 base_height = 9.75
)