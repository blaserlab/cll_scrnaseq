density_umap_pt_timepoint_merged <- 
  bb_var_umap(cds_main[,colData(cds_main)$partition_assignment == "B"], 
              var = "density", 
              facet_by = c("patient","timepoint_merged"), 
              rows = vars(patient), cols = vars(timepoint_merged)) +
  theme(panel.background = element_rect(color = "grey80"))
density_umap_pt_timepoint_merged


