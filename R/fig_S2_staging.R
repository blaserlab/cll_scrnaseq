subpop_umap_pt_timepoint_merged <- 
  bb_var_umap(cds_main[,colData(cds_main)$partition_assignment == "B"], 
              var = "leiden_assignment_binned_renamed", 
              facet_by = c("patient","timepoint_merged"), 
              rows = vars(patient), cols = vars(timepoint_merged)) +
  theme(panel.background = element_rect(color = "grey80"))
subpop_umap_pt_timepoint_merged

