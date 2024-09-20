mono_subpop_umap <-
  bb_var_umap(
    cds_main[, colData(cds_main)$partition_assignment == "Mono"],
    "seurat_l2_leiden_consensus",
    overwrite_labels = T,
    group_label_size = 4,
    foreground_alpha = 0.2,
    rasterize = TRUE
  )

