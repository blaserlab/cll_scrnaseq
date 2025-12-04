
make_fig_S4 <- function(source_all_figs = TRUE) {
  if (source_all_figs) {
    source("R/fig_staging/bcr_sig_hm.R", local = TRUE)
    
  }
  
  save_plot(
    plot_grid(bcr_sig_hm),
    filename = fs::path(network_out, "fig_S4.pdf"),
    base_width = 7.5,
    base_height = 9.75
  )
}
make_fig_S4()