
# subcluster_umap-----------------------------------
umap_subcluster <- 
  bb_var_umap(
    obj = cds_main[,colData(cds_main)$partition_assignment == "B"],
    var = "leiden_assignment_binned_renamed",
    cell_size = 1,
    foreground_alpha = 0.1,
    overwrite_labels = TRUE,
    group_label_size = 4
  ) 
umap_subcluster
