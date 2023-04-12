bb_cellmeta(cds_main) |> glimpse()
bb_var_umap(cds_main, "partition_assignment")
bb_var_umap(cds_main, "leiden_assignment_binned_renamed")
bb_var_umap(cds_main, "seurat_celltype_l1")


bb_var_umap(cds_main, "leiden_l1_assignment")

umap_leiden_l1 <- bb_var_umap(cds_main, "leiden_l1_assignment", value_to_highlight = c("B", "CD4 T", "CD8 T", "NK"), facet_by = c("timepoint_merged_1", "patient_type1"), cols = vars(timepoint_merged_1), rows  = vars(patient_type1), foreground_alpha = 0.2) + panel_border()
save_plot(umap_leiden_l1, filename = fs::path(network_out, "umap_leiden_l1.png"), base_height = 4.0, base_width = 7.5)



