# source("R/dependencies.R")
# source("R/configs.R")
source("R/fig_staging/go_bubbles.R")

fig_S2 <- plot_grid(stressed_goscatter,
                    infl1_goscatter, 
                    infl2_goscatter, 
                    bcr_goscatter, 
                    ncol = 1,
                    align = "v",
                    axis = "lr",
                    labels = c("A", "B", "C", "D")
                    )

save_plot(
 plot = fig_S2,
 filename = "test.png",
 # filename = str_glue("{network_out}/fig_S2.png"),
 base_height = 9.75,
 base_width = 7.5
)
