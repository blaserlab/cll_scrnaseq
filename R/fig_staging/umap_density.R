# B cell umap density faceted--------------------------------
umap_density <- 
  bb_var_umap(
    obj = cds_main,
    # obj = cds_main[,colData(cds_main)$partition_assignment == "B"],
    var = "density",
    sample_equally = TRUE, 
    cell_size = 0.5,
    facet_by = c("patient_type3", "timepoint_merged_1"),
    # rows = vars(patient_type2), 
    # cols = vars(timepoint_merged_2),
    foreground_alpha = 0.2, 
    rasterize = FALSE
  ) +
  theme(panel.background = element_rect(color = "grey80")) +
  theme(legend.justification = "center") +
  labs(color = "Cell\nDensity")
umap_density
