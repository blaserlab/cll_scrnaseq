# source("R/dependencies.R")
# source("R/configs.R")

source_all_figs <- TRUE
# source_all_figs <- FALSE

if (source_all_figs) {
  source("R/fig_staging/cellchat_heatmap.R")
  source("R/fig_staging/heatmap_validation.R")
}

fig_5_top <- plot_grid(
  cellchat_hm,
  labels = c("A")
)

fig_5_bot <- plot_grid(
  cellchat_validation_plots$`MIF | CD8 TEM | CD8 TEM`,
  cellchat_validation_plots$`MIF | B memory | CD14 Mono`,
  ncol = 2,
  rel_widths = c(1,1),
  labels = c("B", "C")
)

fig_5 <- plot_grid(
  fig_5_top,
  fig_5_bot,
  ncol = 1,
  rel_heights = c(1,0.4)
)

save_plot(
  fig_5,
  filename = fs::path(network_out, "fig_5.pdf"),
  base_width = 7.5,
  base_height = 9.75
)
