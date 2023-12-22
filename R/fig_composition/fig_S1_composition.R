# source("R/dependencies.R")
# source("R/configs.R")

source_all_figs <- TRUE
# source_all_figs <- FALSE

if (source_all_figs) {
  source("R/fig_staging/reference_anno_barchart.R")
  source("R/fig_staging/partition_assignment_genebubbles.R")
    
}

fig_S1 <- 
  plot_grid(
    reference_anno_barchart + theme(plot.margin = unit(c(5, 5, 0, 0), "mm")),
    partition_assignment_genebubbles,
    nrow = 2,
    align = "v",
    axis = "l",
    labels = c("A", "B"),
    rel_heights = c(1,2))

save_plot(
 plot = fig_S1,
 # filename = "fig_s1.png",
 filename = fs::path(network_out, "fig_S1.pdf"),
 base_width = 7.5,
 base_height = 9.75
)