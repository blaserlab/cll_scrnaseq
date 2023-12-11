# source("R/dependencies.R")
# source("R/configs.R")
source("R/fig_staging/pt_characteristics_hm.R")

fig_S2 <- plot_grid(pt_characteristics_hm)

save_plot(
 plot = fig_S2,
 filename = "test.png",
 # filename = str_glue("{network_out}/fig_S2.png"),
 base_height = 7.5,
 base_width = 9.75
)
