colData(cds_main)$patient <- factor(colData(cds_main)$patient, levels = str_sort(unique(colData(cds_main)$patient),numeric = T))
colData(cds_main)$timepoint_merged <- recode(colData(cds_main)$timepoint, 
                                             "3yrs" = "3yrs|btk_clone",
                                             "btk_clone" = "3yrs|btk_clone",
                                             "5yrs" = "5yrs|relapse",
                                             "relapse" = "5yrs|relapse",
                                             )

density_umap_pt_timepoint_merged <- 
  bb_var_umap(cds_main[,colData(cds_main)$partition_assignment_2 == "B"], 
              var = "density", 
              facet_by = c("patient","timepoint_merged"), 
              rows = vars(patient), cols = vars(timepoint_merged)) +
  theme(panel.background = element_rect(color = "grey80"))
