# global umap with partition labels ----------------------------------
umap_partitions <- bb_var_umap(
  cds = cds_main,
  var = "partition_assignment_2",
  overwrite_labels = T,
  group_label_size = 4,
  foreground_alpha = 0.1
)

umap_bcell_leiden <- 
  bb_var_umap(cds_main[, colData(cds_main)$leiden_assignment_1 %in% c("BTK cluster", "MRD1 cluster", "MRD2 cluster")],
              var = "leiden", 
              overwrite_labels = T,
              foreground_alpha = 0.1)

umap_binned_leiden <-
  bb_var_umap(cds_main[, colData(cds_main)$leiden_assignment_1 %in% c("BTK cluster", "MRD1 cluster", "MRD2 cluster")],
              var = "leiden_assignment_1", 
              overwrite_labels = T,
              foreground_alpha = 0.1)