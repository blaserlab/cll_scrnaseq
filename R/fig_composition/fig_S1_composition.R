source("R/dependencies.R")
source("R/configs.R")
source("R/fig_staging/reference_anno_barchart.R")
source("R/fig_staging/partition_assignment_genebubbles.R")


fig_S1_top <- 
  plot_grid(
    reference_anno_barchart,
    ncol = 1,
    labels = c("A")
  )


fig_S1_bot <- 
  plot_grid(partition_assignment_genebubbles,
            labels = c("B"))



fig_S1 <- 
  plot_grid(
    fig_S1_top,
    fig_S1_bot,
    nrow = 2,
    rel_heights = c(1,2))

save_plot(
 plot = fig_S1,
 filename = "test.png",
 # filename = str_glue("{network_out}/fig_S1.png"),
 base_width = 7.5,
 base_height = 9.75
)