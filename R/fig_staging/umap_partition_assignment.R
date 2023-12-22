
umap_partition_assignment <- bb_var_umap(
  obj = cds_main,
  var = "partition_assignment",
  foreground_alpha = 0.05,
  rasterize = TRUE,
  overwrite_labels = TRUE,
  palette = experimental_group_palette
)
umap_partition_assignment
