# source("R/dependencies.R")
# source("R/configs.R")

source_all_figs <- TRUE
# source_all_figs <- FALSE

if (source_all_figs) {
  source("R/fig_staging/go_bubbles.R")
    
}


fig_S3 <- plot_grid(stressed_goscatter,
                    infl1_goscatter, 
                    infl2_goscatter, 
                    bcr_goscatter, 
                    ncol = 1,
                    align = "v",
                    axis = "lr",
                    labels = c("A", "B", "C", "D")
                    )

save_plot(
 plot = fig_S3,
 # filename = "fig_s3.png",
 filename = fs::path(network_out, "fig_S3.pdf"),
 base_height = 9.75,
 base_width = 7.5
)
