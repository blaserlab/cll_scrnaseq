# source("R/dependencies.R")
# source("R/configs.R")

source_all_figs <- TRUE
# source_all_figs <- FALSE

if (source_all_figs) {
  source("R/fig_staging/pt_characteristics_hm.R")
  source("R/fig_staging/cll_flow.R")
    
}

fig_S2 <- plot_grid(
  pt_characteristics_hm,
  plot_grid(cll_flow_fig, NULL, ncol = 2, labels = c("B", "")),
  ncol = 1,
  rel_heights = c(2,1), labels = c("A", ""))

save_plot(
 plot = fig_S2,
 filename = fs::path(network_out, "fig_S2.pdf"),
 base_width = 7.5,
 base_height = 9.75
)
