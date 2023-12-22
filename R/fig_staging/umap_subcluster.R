
# subcluster_umap-----------------------------------
umap_subcluster <- 
  bb_var_umap(
    obj = cds_main[,colData(cds_main)$partition_assignment == "B"],
    var = "leiden_comparison_renamed",
    cell_size = 1,
    foreground_alpha = 0.1,
    overwrite_labels = TRUE,
    palette = experimental_group_palette,
    rasterize = TRUE
  ) 
umap_subcluster
