source("R/configs.R")
source("R/fig_1_staging.R")

fig_1_top <-
  plot_grid(
    umap_density,
    NULL,
    ncol = 2,
    rel_widths = c(4,1),
    labels = c("A","")
    )
    
fig_1_mid <- 
  plot_grid(
    volcano_BTK_bcells,
    volcano_MRD1_BTK,
    volcano_MRD2_BTK,
    labels = c("B","C","D"),
    align = "h",
    axis = "b",
    ncol = 3,
    rel_widths = c(1,1,1)
  )



fig_1_bottom <- 
  plot_grid(
    module_heatmap_bcells,
    subpop_gene_dotplot,
    ncol = 2,
    rel_widths = c(1,1),
    labels = c("E","F")
  )

fig_1 <- 
  plot_grid(
    fig_1_top,
    fig_1_mid,
    fig_1_bottom,
    nrow = 3,
    rel_heights = c(1.2,1,1.5)
  )

save_plot(
 plot = fig_1,
 filename = str_glue("{network_out}/fig_1.png"),
 base_width = 7.5,
 base_height = 9.75
)