source("00_packages_functions.R")

pt_timepoint_vec <- unique(colData(cds_aligned)$pt_timepoint)

pt_timepoint_facet_plot <- custom_variable_plot(
  cds = cds_aligned,
  var = "pt_timepoint",
  foreground_alpha = 0.4,
  value_to_highlight = pt_timepoint_vec,
  palette = rep("red",times = 12),
  legend_pos = "none"
) + theme(panel.background = element_rect(color = "grey80"))+
  facet_grid(cols = vars(timepoint), rows = vars(pt)) +
  theme(strip.background = element_blank())
save_plot(
  pt_timepoint_facet_plot,
  filename = "plots_out/pt_timepoint_facet_plot.pdf",
  base_height = 6.5,
  base_width = 7
)

# make a folder
if (!dir.exists("plots_out/single_pt_timepoint_plots")) {
  dir.create("plots_out/single_pt_timepoint_plots")
}

#singles
pt_timepoint_plots<-lapply(
  X = seq_along(pt_timepoint_vec),
  FUN = custom_variable_plot,
  cds = cds_aligned,
  var = "pt_timepoint", 
  value_to_highlight = pt_timepoint_vec, 
  foreground_alpha = 0.4, 
  legend_pos = "none", 
  cell_size = 0.5, 
  legend_title = NULL, 
  plot_title = pt_timepoint_vec, 
  outfile = paste0("single_pt_timepoint_plots/",pt_timepoint_vec), 
  h = 4, 
  w = 4.4, 
  palette = "red"
)


save.image.pigz("cll_scrnaseq.RData",n.cores = 39)
