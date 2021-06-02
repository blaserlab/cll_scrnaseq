source("R/configs.R")
source("R/fig_S1_staging.R")

fig_S1 <- 
  plot_grid(
    umap_partitions,NULL,
    umap_bcell_leiden,NULL,
    umap_binned_leiden,NULL,
    nrow = 3,
    rel_widths = c(1.5,1),
    labels = c("A","","B","","",""))

save_plot(
 plot = fig_S1,
 filename = str_glue("{network_out}/fig_S1.png"),
 base_width = 7.5,
 base_height = 9.75
)