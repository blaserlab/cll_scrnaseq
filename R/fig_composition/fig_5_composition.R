
make_fig_5 <- function(source_all_figs = TRUE) {
  if (source_all_figs) {
    source("R/fig_staging/cellchat_figs.R", local = TRUE)
    source("R/fig_staging/cellchat_flow.R", local = TRUE)
  }
  
  fig_5_left <- plot_grid(cellchat_hm,  labels = c("A"))
  
  fig_5_right <- plot_grid(
    cellchat_validation_plots$MIF_CD74_CD44$plot,
    cellchat_flow_figs$MIF$plot,
    cellchat_validation_plots$LGALS9_CD45$plot,
    cellchat_flow_figs$LGALS$plot,
    ncol = 1,
    align = "v", 
    axis = "l",
    labels = c("B", "C", "D", "E")
  )
  
  fig_5 <- plot_grid(fig_5_left,
                     fig_5_right,
                     ncol = 2,
                     rel_widths = c(2, 1))
  
  save_plot(
    fig_5,
    filename = fs::path(network_out, "fig_5.pdf"),
    base_width = 7.5,
    base_height = 9.75
  )
  
}
make_fig_5()