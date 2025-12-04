

make_fig_2 <- function(source_all_figs = TRUE) {
  if (source_all_figs) {
    source("R/fig_staging/module_expression_hm.R", local = TRUE)
    source("R/fig_staging/mod4_enrichment_plot.R", local = TRUE)
    source("R/fig_staging/mod4_agg_expr.R", local = TRUE)
    source("R/fig_staging/mod4_genebub.R", local = TRUE)
  }
  
  module_heatmap_bcells_alt <-
    plot_grid(module_heatmap_bcells) 
  
  fig_2_top <-
    plot_grid(module_heatmap_bcells_alt,mod4_agg_expr_violin,
              ncol = 2,
              labels = c("A", "B"))
  
  fig_2_mid <-
    plot_grid(
      mod4_genebub,
      ncol = 1,
      labels = c("D")
    )
  
  
  fig_2_bottom <-
    plot_grid(mod4_enrichment_plot, labels = "C")
  
  fig_2 <- plot_grid(
    fig_2_top,
    fig_2_bottom,
    fig_2_mid,
    ncol = 1,
    rel_heights = c(0.5, 0.6, 0.5)
  )
  
  save_plot(
    plot = fig_2,
    # filename = "fig_2.png",
    filename = fs::path(network_out, "fig_2.pdf"),
    base_width = 7.5,
    base_height = 9.75
  )
}
make_fig_2()