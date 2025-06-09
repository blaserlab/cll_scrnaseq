# source("R/dependencies.R")
# source("R/configs.R")

source_all_figs <- TRUE
# source_all_figs <- FALSE

if (source_all_figs) {
  source("R/fig_staging/bcr_sig_hm.R")
  
}



save_plot(
  plot_grid(bcr_sig_hm),
  filename = fs::path(network_out, "fig_S2.5.pdf"),
  base_width = 7.5,
  base_height = 9.75
)
