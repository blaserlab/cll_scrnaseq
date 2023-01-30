fig_3_left <- 
  align_plots(
    tcell_subpop_gene_umap,
    tcell_subpop_umap,
    major_tcell_cluster_umap,
    align = "v", 
    axis = "l"
  )


fig_3_top <-
  plot_grid(fig_3_left[[1]],
            ncol = 1,
            rel_widths = c(1.0),
            labels = "A")

fig_3_2 <-
  plot_grid(
    fig_3_left[[2]],
    tcell_subpop_gene_dotplot,
    labels = c("B", "C"),
    ncol = 2, 
    rel_widths = c(0.8, 1.0),
    align = "h",
    axis = "b"
  )

fig_3_3 <-
  plot_grid(
    fig_3_left[[3]],
    treg_pct_plot,
    labels = c("D", "E"),
    align =  "h",
    axis = "b",
    ncol = 2,
    rel_widths = c(0.8,1.0)
  )

fig_3 <- 
  plot_grid(
    fig_3_top,
    fig_3_2,
    fig_3_3,
    align = "v",
    axis = "l",
    nrow = 3, 
    rel_heights = c(1.5, 1.0, 1.0)
  )

save_plot(fig_3, filename = str_glue("{network_out}/fig_3.png"), base_width = 7.5, base_height = 9.75)