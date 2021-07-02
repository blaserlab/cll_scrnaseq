source("R/configs.R")
source("R/fig_S2_staging.R")

fig_S2 <- density_umap_pt_timepoint_merged 

save_plot(
 plot = fig_S2,
 filename = str_glue("{network_out}/fig_S2.png"),
 base_width = 7.5,
 base_height = 9.75
)
