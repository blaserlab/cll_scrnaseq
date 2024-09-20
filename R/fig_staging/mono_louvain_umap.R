mono_louvain_umap <-
  bb_var_umap(
    cds_main[, colData(cds_main)$partition_assignment %in% c("Mono")],
    "louvain",
    overwrite_labels = T,
    group_label_size = 4,
    foreground_alpha = 0.2,
    rasterize = TRUE
  )
