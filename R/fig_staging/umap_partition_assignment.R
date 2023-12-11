
umap_partition_assignment <- bb_var_umap(
  obj = cds_main,
  var = "partition_assignment",
  foreground_alpha = 0.05,
  rasterize = TRUE,
  overwrite_labels = TRUE
)
umap_partition_assignment + bb_var_umap(cds_main, "seurat_celltype_l1", rasterize = TRUE)
