tnk_leiden_umap <-
  bb_var_umap(
    cds_main[, colData(cds_main)$partition_assignment %in% c("T", "NK")],
    "leiden",
    overwrite_labels = T,
    group_label_size = 4,
    foreground_alpha = 0.2,
    rasterize = TRUE
  )
