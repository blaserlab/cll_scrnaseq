fig_2_top <- plot_grid(
  tcell_subpop_umap,
  tcell_subpop_genebub,
  ncol = 2,
  align = "h", 
  axis = "b",
  rel_widths = c(1,1)
)

fig_2_mid <- plot_grid(
  cluster_enrichment_barchart,
  treg_pct_plot,
  tcr_diversity_plot,
  ncol = 3,
  rel_widths = c(1,1,1)
)

fig_2_bot <- plot_grid(
  exh_genebub,
  NULL,
  ncol =  2,
  rel_widths = c(1,1)
)

fig_2 <- plot_grid(
  fig_2_top,
  fig_2_mid,
  fig_2_bot,
  ncol = 1,
  rel_heights = c(1,1,1)
)

save_plot(
  fig_2,
  filename = fs::path(network_out, "fig_2.png"),
  base_width = 7.5,
  base_height = 9.75
)
