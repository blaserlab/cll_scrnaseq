# source("R/dependencies.R")
# source("R/configs.R")

source_all_figs <- TRUE
# source_all_figs <- FALSE

if (source_all_figs) {
  source("R/fig_staging/module_expression_hm.R")
  source("R/fig_staging/mod4_enrichment_plot.R")
  source("R/fig_staging/mod4_agg_expr.R")
}

module_heatmap_bcells_alt <-
  plot_grid(module_heatmap_bcells) + theme(plot.margin = unit(c(0.25, 1.5, 0.25, 0.5), "in"))

fig_2_top <-
  plot_grid(module_heatmap_bcells_alt,
            ncol = 1,
            labels = c("A"))

fig_2_mid <-
  plot_grid(
    mod4_agg_expr_violin,
    NULL,
    ncol = 2,
    rel_widths = c(3, 1),
    labels = c("B")
  )


fig_2_bottom <-
  plot_grid(mod4_enrichment_plot,
            labels = "C")

fig_2 <- plot_grid(fig_2_top,
                   fig_2_mid,
                   fig_2_bottom,
                   ncol = 1,
                   rel_heights = c(0.6, 0.4, 0.6))

save_plot(
  plot = fig_2,
  # filename = "fig_2.png",
  filename = fs::path(network_out, "fig_2.pdf"),
  base_width = 7.5,
  base_height = 9.75
)
