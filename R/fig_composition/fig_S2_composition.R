# source("R/dependencies.R")
# source("R/configs.R")

source_all_figs <- TRUE
# source_all_figs <- FALSE

if (source_all_figs) {
  source("R/fig_staging/pt_characteristics_hm.R")
    
}

fig_S2 <- plot_grid(pt_characteristics_hm)

save_plot(
 plot = fig_S2,
 # filename = "fig_s2.png",
 filename = fs::path(network_out, "fig_S2.pdf"),
 base_height = 7.5,
 base_width = 9.75
)
