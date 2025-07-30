# source("R/dependencies.R")
# source("R/configs.R")

source_all_figs <- TRUE
# source_all_figs <- FALSE

if (source_all_figs) {
  source("R/fig_staging/tsubtype_genebub.R")
  
  
}



save_plot(
  plot_grid(tsubtype_genebub),
  filename = fs::path(network_out, "fig_S6.pdf"),
  base_width = 7.5,
  base_height = 9.75
)