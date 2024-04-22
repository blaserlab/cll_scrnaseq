tnk_seurat_celltypes <-
  bb_var_umap(
    cds_main[, colData(cds_main)$partition_assignment %in% c("T", "NK")],
    "seurat_celltype_l2",
    overwrite_labels = T,
    group_label_size = 4,
    foreground_alpha = 0.2,
    rasterize = TRUE
  )
tnk_seurat_celltypes
